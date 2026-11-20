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
//! \brief Closed-form (non-GJK/EPA) contact-detection algorithms from published papers, for the shape pairs
//!        SolidCylinderContactDetection (ParticleKernels.h) otherwise resolves via the generic GJK/EPA fallback:
//!        Box x Box and Cylinder x Cylinder.
//!
//! Thin aggregator of the four themed headers:
//!   - ClosedFormContactsInfra.h: the shared ContactResult return type.
//!   - ClosedFormContactsBoxPoint.h: Box x Box (Gottschalk/Lin/Manocha 1996 SAT + this codebase's own
//!     contact-point construction).
//!   - ClosedFormContactsCylinderPlaneGeometry.h: CylinderState, the quartic/ellipse/rim-crossing
//!     solvers, Cylinder x Plane, Cylinder x Sphere.
//!   - ClosedFormContactsCylinderCylinder.h: Cylinder x Cylinder (depends on the geometry file above).
//!
//! Use in the simulation is opt-in: SolidCylinderContactDetection (ParticleKernels.h, used by ComputeContacts) calls
//! boxBoxSAT / cylinderCylinderContact / cylinderCylinderResultIsUncertain only if
//! ModelParameters.useAnalyticalContactDetection is set (default false, GJK/EPA only; the diagnostic
//! shadowAnalyticCylinderLog mode evaluates cylinderCylinderContact alongside GJK/EPA without using it). A result
//! flagged by cylinderCylinderResultIsUncertain still falls back to GJK/EPA unless
//! ModelParameters.usePureAnalyticalCylinder is also set, which trusts the flagged analytical answer with no safety
//! net; see SolidCylinderContactDetection's Cylinder x Cylinder operator() for the exact tiering (including the
//! opt-in refineCylinderContacts path). PaperAnalyticalContactTest.cpp compares
//! these closed forms against SolidCylinderContactDetection's GJK/EPA results standalone.
//!
//! -------------------------------------------------------------------------------------------------
//! Box x Box: Gottschalk, S., Lin, M. C., & Manocha, D. (1996). "OBBTree: A Hierarchical Structure
//! for Rapid Interference Detection." SIGGRAPH 1996 Proceedings, 171-180. Section 5, "Fast Overlap
//! Test for OBBs" (the 15-axis separating-axis-theorem test).
//!
//! The paper derives only a BOOLEAN overlap test. The minimum-translation-vector extension used here for a
//! penetration depth and normal (track the axis with the SMALLEST overlap instead of exiting on the first
//! separating axis) is the standard textbook extension of the same theorem (e.g. Ericson, "Real-Time Collision
//! Detection", 2004, ch. 4), not a result of the 1996 paper. The contact-POINT construction is not from the paper
//! either. A flush face-face contact has no unique contact point: ContactDetectionTest.cpp's testBoxBoxContact therefore
//! asserts only the normal-axis component of the GJK/EPA point (pybullet picks a different, equally valid point on the
//! same shared face).
//!
//! -------------------------------------------------------------------------------------------------
//! Cylinder x Cylinder / Cylinder x Plane: Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., &
//! Wassgren, C. (2010). "Cylindrical object contact detection for use in discrete element method
//! simulations. Part I - Contact detection algorithms." Chemical Engineering Science, 65(22),
//! 5852-5862. All "(eq. N)" comments below refer to this paper's equation numbers.
//!
//! Validation caveat (companion Part II paper, Chemical Engineering Science 65(22), 5863-5871): Part II validates
//! the cylinder-PLANE algorithm directly against a closed-form rebound formula (single cylinder bouncing off a
//! wall). It does NOT validate cylinder-CYLINDER contact geometry (overlap/normal/point) against direct ground
//! truth; its multi-particle experiments (fill height and bulk solid fraction of 250 cylinders dropped into a
//! container, dynamic angle of repose in a rotating drum) compare only AGGREGATE bulk measurements. "This is the
//! published algorithm" is therefore a stronger
//! claim than "its cylinder-cylinder branch has been independently ground-truthed".
//!
//! Corrections: Guo, Y., Wassgren, C., Ketterhagen, W., Hancock, B., & Curtis, J. (2012). "Some
//! computational considerations associated with discrete element modeling of cylindrical particles."
//! Powder Technology, 228, 193-198. A follow-up by an overlapping author team documenting classification bugs in
//! Kodam et al.'s face-edge/band-edge/band-band-skewed dispatch ("undetected or wrongly detected contacts can lead
//! to large overlaps") and a scenario Kodam et al. never defined (their Fig. 8 "special contact").
//! ClosedFormContactsCylinderCylinder.h implements Guo et al.'s corrected criteria (see rimVsFaceOrBand() and the
//! band-band-skewed exclusion) in place of Kodam et al.'s originals; Kodam's original dispatch gives false negatives
//! for overlapping, obliquely tilted pairs that GJK/EPA detects (PaperAnalyticalContactTest.cpp).
//!
//! Two sub-problems (band-edge, eq. 74-77; edge-edge point location, eq. 82-85) reduce to a fourth-order
//! polynomial in sin(theta). closestThetaOnCircleToLine (band-edge) and findRimCrossing (edge-edge point A) instead
//! derive quartics in t = tan(theta/2) directly from the geometry (algebra in their own comments) and solve them
//! with Ferrari's method (solveQuarticFerrari, O(1), no iteration), the default; QuarticSolverComparisonTest.cpp
//! found it as accurate as Durand-Kerner iteration at a fraction of the cost. Durand-Kerner and Aberth-Ehrlich remain for comparison (see QuarticSolverMethod).
//! closestThetaOnCircleToPoint (findRimCrossing's optional fast reference angle) has a simpler exact answer, a single
//! atan2, since its objective is a pure first-harmonic trig function of theta. Both closed forms match a
//! golden-section reference minimisation to floating-point precision across thousands of randomised
//! configurations.
//======================================================================================================================

#pragma once

#include "ClosedFormContactsInfra.h"
#include "ClosedFormContactsBoxPoint.h"
#include "ClosedFormContactsCylinderPlaneGeometry.h"
#include "ClosedFormContactsCylinderCylinder.h"
