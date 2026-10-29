//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file ClosedFormContacts.h
//! \brief Closed-form (non-GJK/EPA) contact-detection algorithms taken directly from published
//!        papers, for the shape pairs SolidCylinderContactDetection (ParticleKernels.h) currently
//!        resolves via the generic GJK/EPA fallback: Box x Box and Cylinder x Cylinder.
//!
//! [Claude, 2026-10-29] Renamed from PaperAnalyticalContacts.h (namespace paper_analytical ->
//! closed_form_contacts) and split into 4 files by theme -- this header is now a thin aggregator,
//! kept so every external #include "PaperAnalyticalContacts.h" site needed only a filename update,
//! not a restructure:
//!   - ClosedFormContactsInfra.h: the shared ContactResult return type.
//!   - ClosedFormContactsBoxPoint.h: Box x Box (Gottschalk/Lin/Manocha 1996 SAT + this codebase's own
//!     contact-point construction and vertex-only multi-point contact).
//!   - ClosedFormContactsCylinderPlaneGeometry.h: CylinderState, the quartic/ellipse/rim-crossing
//!     solvers, Cylinder x Plane, Cylinder x Sphere.
//!   - ClosedFormContactsCylinderCylinder.h: Cylinder x Cylinder (depends on the geometry file above).
//!
//! USED BY THE SIMULATION, opt-in: ParticleKernels.h's ComputeContacts includes this header and calls
//! cylinderCylinderContact/cylinderCylinderResultIsUncertain directly, gated behind the runtime prm
//! flag ModelParameters.useAnalyticalContactDetection (default false -- GJK/EPA-only unless a scenario
//! explicitly turns this on). When on, a result cylinderCylinderResultIsUncertain flags still falls
//! back to GJK/EPA UNLESS ModelParameters.usePureAnalyticalCylinder is ALSO set, which trusts the flagged
//! analytical answer with no safety net -- see ComputeContacts::operator()'s own comment for the exact
//! tiering. PaperAnalyticalContactTest.cpp additionally compares these closed forms against
//! SolidCylinderContactDetection's GJK/EPA results directly, standalone from the production dispatch.
//!
//! -------------------------------------------------------------------------------------------------
//! Box x Box: Gottschalk, S., Lin, M. C., & Manocha, D. (1996). "OBBTree: A Hierarchical Structure
//! for Rapid Interference Detection." SIGGRAPH 1996 Proceedings, 171-180. Section 5, "Fast Overlap
//! Test for OBBs" (the 15-axis separating-axis-theorem test).
//!
//! The paper itself only derives a BOOLEAN overlap test (is there a separating axis, yes/no). The
//! minimum-translation-vector extension used here to also produce a penetration depth and contact
//! normal -- track the axis with the SMALLEST overlap instead of exiting on the first separating
//! axis found -- is the standard textbook extension of the same separating-axis theorem (see e.g.
//! Ericson, "Real-Time Collision Detection", 2004, ch. 4), not an additional result the 1996 paper
//! itself derives. A supplementary contact-POINT construction (also not from the paper) is provided
//! too, restricted to the one well-defined axis component -- the same restriction
//! ContactDetectionTest.cpp's testBoxBoxContact already applies to the GJK/EPA result, since a
//! flush face-face contact has no unique contact point in the first place (confirmed there by
//! cross-checking against pybullet, which picks a different, equally valid, point on the same
//! shared face).
//!
//! -------------------------------------------------------------------------------------------------
//! Cylinder x Cylinder / Cylinder x Plane: Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., &
//! Wassgren, C. (2010). "Cylindrical object contact detection for use in discrete element method
//! simulations. Part I - Contact detection algorithms." Chemical Engineering Science, 65(22),
//! 5852-5862. All "(eq. N)" comments below refer to this paper's equation numbers.
//!
//! Validation caveat (from the companion Part II paper, Chemical Engineering Science 65(22),
//! 5863-5871): Part II validates the cylinder-PLANE algorithm directly against a closed-form
//! analytical rebound formula (single cylinder bouncing off a wall). It does NOT validate
//! cylinder-CYLINDER contact geometry (overlap/normal/point) against any direct ground truth --
//! the only cylinder-cylinder-contact-bearing experiment (250 cylinders packed in a container) is
//! validated against an AGGREGATE bulk measurement (fill height, bulk solid fraction), not against
//! individual contact geometry. So "this is the published algorithm" is a stronger claim than "this
//! algorithm's cylinder-cylinder branch has been independently ground-truthed" -- keep that
//! distinction in mind when reading the test comparisons below.
//!
//! CORRECTIONS: Guo, Y., Wassgren, C., Ketterhagen, W., Hancock, B., & Curtis, J. (2012). "Some
//! computational considerations associated with discrete element modeling of cylindrical particles."
//! Powder Technology, 228, 193-198. An overlapping author team's own follow-up to Kodam et al. (2010),
//! documenting real classification bugs in the original face-edge/band-edge/band-band-skewed dispatch
//! logic ("undetected or wrongly detected contacts can lead to large overlaps") and one whole scenario
//! Kodam et al. never defined at all (their Fig. 8 "special contact"). This file now implements Guo et
//! al.'s corrected criteria (see rimVsFaceOrBand() and the band-band-skewed exclusion below) in place
//! of Kodam et al.'s original ones.
//!
//! HISTORY (kept for context, since this is exactly the kind of thing worth being able to point to):
//! an earlier version of this file used Kodam et al.'s ORIGINAL (uncorrected) face-edge/band-edge
//! dispatch logic and was demonstrated by PaperAnalyticalContactTest.cpp to produce false negatives
//! -- reporting no contact for genuinely overlapping, oblique-tilted cylinder pairs that GJK/EPA
//! correctly detects. Applying Guo et al.'s corrections (this version) targets exactly that class of
//! scenario-boundary misclassification; see PaperAnalyticalContactTest.cpp's own comments for which
//! specific previously-failing cases now pass and which, if any, still don't.
//!
//! Two sub-problems in the paper (band-edge, eq. 74-77; edge-edge point location, eq. 82-85) reduce
//! to a fourth-order polynomial in sin(theta). An earlier version of this file avoided hand-
//! transcribing that quartic (risk of an undetected sign/coefficient error with no compiler-
//! verified reference) and instead solved the equivalent 1D minimization numerically: a coarse
//! uniform scan over theta followed by golden-section refinement (up to ~1500 evaluations per
//! candidate contact). That numeric approach turned out to be a genuine performance problem at
//! scale -- dense scenarios with many simultaneous Cylinder x Cylinder contacts ran ~6x slower per
//! step than an equivalent Box x Box scenario purely because of this search.
//!
//! This version replaces it with the actual closed form: closestThetaOnCircleToLine derives the
//! quartic coefficients directly from the geometry (see that function's own comment for the
//! algebra) and solves them via the Durand-Kerner method (solveQuarticDurandKerner) -- chosen over
//! Ferrari's method specifically because it has no case-branching on a resolvent cubic's
//! discriminant to get subtly wrong; every root, real or complex, is refined by the same fixed
//! iteration loop. closestThetaOnCircleToPoint (used by the edge-edge alternating-projection loop)
//! turned out to have an even simpler exact answer -- a single atan2, no polynomial at all, since
//! its objective is a pure first-harmonic trig function of theta. Both closed forms were validated
//! against the previous golden-section implementation across thousands of randomized
//! configurations before replacing it (matching to floating-point precision in every trial), the
//! "compiler-verified reference" the original numeric approach didn't have.
//======================================================================================================================

#pragma once

#include "ClosedFormContactsInfra.h"
#include "ClosedFormContactsBoxPoint.h"
#include "ClosedFormContactsCylinderPlaneGeometry.h"
#include "ClosedFormContactsCylinderCylinder.h"
