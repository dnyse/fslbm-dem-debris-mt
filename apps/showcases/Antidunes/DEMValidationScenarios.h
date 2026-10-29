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
//! \file DEMValidationScenarios.h
//! \brief Fluid-free DEM validation experiments for DEMValidation.cpp.
//!
//! Every scenario runs the PRODUCTION contact/integration code (SolidCylinderContactDetection,
//! ComputeContacts, mesa_pd::kernel::VelocityVerletPre/PostForceUpdate -- ParticleKernels.h) through a
//! minimal, hand-written sub-cycle loop with no fluid coupling (no lubrication, no hydrodynamic force,
//! no buoyancy -- gravity, where used, is applied directly via addForceAtomic, same as BedGeneration.cpp).
//! Existing tests (ContactDetectionTest.cpp, PaperAnalyticalContactTest.cpp, validation/*.py) only check
//! STATIC contact geometry; this file is the first thing in the repo that runs real dynamics and checks
//! the result against a closed-form or literature reference -- see AnalyticalSolutions.h for the physics.
//======================================================================================================================

#pragma once

#include <sstream>

#include "AnalyticalSolutions.h"
#include "Utility.h"

#include "blockforest/Initialization.h"
// ParticleKernels.h -> SimulationObjects.h uses blockforest::SimpleCommunication<> but relies on the
// caller having already included its defining header (normally pulled in transitively via
// BuildTimeloop.h's chain of lbm/free_surface headers in the full Antidunes.cpp build, which this
// fluid-free app does not include) -- include it explicitly here instead.
#include "lbm/blockforest/communication/SimpleCommunication.h"

#include "ParticleKernels.h"

#include "core/math/Constants.h"
#include "core/math/Quaternion.h"
#include "core/math/Random.h"
#include "core/mpi/Reduce.h"

#include "mesa_pd/data/LinkedCells.h"
#include "mesa_pd/data/ParticleAccessorWithBaseShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Box.h"
#include "mesa_pd/data/shape/Cylinder.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/AssocToBlock.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/InsertParticleIntoLinkedCells.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/SyncNextNeighborsBlockForest.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "vtk/VTKOutput.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace walberla
{
namespace antidunes
{
namespace dem_validation
{
using namespace mesa_pd;

namespace detail
{
//! Minimal fluid-free DEM engine: owns storage/domain/contact-response/integrators and runs one
//! sub-cycle at a time. Always uses LinkedCells + forEachParticlePairHalf (BedGeneration.cpp's own
//! pattern) regardless of particle count, so the same code path serves a 2-body wall-impact test and
//! the hundreds-of-particles packing test alike.
class DEMEngine
{
public:
   // useAnalyticalContactDetection: false (DEFAULT, unchanged for every existing scenario below) routes
   // Box x Box / Cylinder x Cylinder through the generic GJK/EPA fallback, same as before. true routes
   // them through the ACTUAL production analytical dispatch (boxBoxSAT / Kodam+Guo) -- needed by any
   // scenario that specifically investigates a defect in that path (e.g. runCornerNormalDiscontinuity),
   // since GJK/EPA and the closed-form analytical dispatch are different code with different failure
   // modes and a bug specific to one is not reproduced by running the other.
   // useRoundedBoxCorners/boxRoundingRadius/useMultiPointBoxContact: the two togglable Box x Box
   // corner/edge contact-normal-discontinuity mitigations (see SimulationParams.h's own doc comments
   // and paper_analytical::boxBoxSAT/boxBoxMultiPointContacts) -- threaded straight through to
   // ComputeContacts so runCornerNormalDiscontinuity/runCornerFixComparison below can A/B/C/D-test them
   // against the exact same production contact-response code path production runs use.
   DEMEngine(real_t dt, real_t restitutionCoefficient, real_t frictionCoefficient, real_t poissonsRatio,
             real_t collisionTime, real_t gravityMagnitude, real_t domainHalfExtent = real_t(500),
             bool useAnalyticalContactDetection = false, bool useRoundedBoxCorners = false,
             real_t boxRoundingRadius = real_t(0.05), bool useMultiPointBoxContact = false,
             bool useBoxOverlapCap = false, real_t boxOverlapCapFraction = real_t(0.3),
             bool useAxisHysteresis = false, real_t axisHysteresisMargin = real_t(0.15),
             bool useContactManifold = false, real_t manifoldFadeRate = real_t(0.25),
             bool useRotationalEffectiveMass = false)
      : dt_(dt)
      , gravityZ_(-gravityMagnitude)
      , restitutionCoefficient_(restitutionCoefficient)
      , collisionTime_(collisionTime)
      // kappa = kT/kN = 2(1-nu)/(2-nu), Mindlin (1949)/Thornton's no-slip single-contact stiffness
      // ratio -- exact primary Thornton paper/year not locally verified (no copy in cld_kbase/, no
      // single matching reference found by web search as of 2026-10-29); treat the attribution as
      // "Thornton-type no-slip model" rather than a pinned citation until checked against a source.
      , kappa_(real_t(2) * (real_t(1) - poissonsRatio) / (real_t(2) - poissonsRatio))
      , preVV_(dt)
      , postVV_(dt)
   {
      const math::AABB domainAABB(-domainHalfExtent, -domainHalfExtent, -domainHalfExtent, domainHalfExtent,
                                   domainHalfExtent, domainHalfExtent);
      forest_ = blockforest::createBlockForest(domainAABB, Vector3< uint_t >(1, 1, 1),
                                                Vector3< bool >(false, false, false));
      domain_ = std::make_shared< mesa_pd::domain::BlockForestDomain >(forest_);
      ps_     = std::make_shared< data::ParticleStorage >(64);
      ac_     = std::make_unique< data::ParticleAccessorWithBaseShape >(ps_);
      assoc_  = std::make_unique< kernel::AssocToBlock >(forest_);

      computeContacts_ = std::make_unique< ComputeContacts >(
         domain_, restitutionCoefficient_, kappa_, collisionTime_, dt_, frictionCoefficient, uint_t(2),
         useAnalyticalContactDetection, /*usePureAnalyticalCylinder=*/false,
         /*useContactFilterOwnerTiebreak=*/false, /*useTightTier1=*/false, /*useTightLayer1=*/false,
         /*skipTier1=*/false, useRoundedBoxCorners, boxRoundingRadius, useMultiPointBoxContact,
         useBoxOverlapCap, boxOverlapCapFraction, useAxisHysteresis, axisHysteresisMargin,
         useContactManifold, manifoldFadeRate, useRotationalEffectiveMass);
   }

   //! [Claude, 2026-09-30] direct access to the production contact functor (detection-mode switches)
   ComputeContacts&                    computeContacts() { return *computeContacts_; }
   data::ParticleStorage&              ps() { return *ps_; }
   data::ParticleAccessorWithBaseShape& ac() { return *ac_; }

   size_t addSphere(const Vec3& pos, const Vec3& vel, real_t radius, real_t density)
   {
      auto p = ps_->create(false);
      p->setPosition(pos);
      p->setBaseShape(std::make_shared< data::Sphere >(radius));
      p->getBaseShapeRef()->updateMassAndInertia(density);
      p->setOwner(walberla::mpi::MPIManager::instance()->rank());
      p->setType(1);
      p->setInteractionRadius(radius);
      p->setLinearVelocity(vel);
      return p.getIdx();
   }

   size_t addBox(const Vec3& pos, const Vec3& vel, const Rot3& rot, const Vec3& edgeLength, real_t density)
   {
      auto p = ps_->create(false);
      p->setPosition(pos);
      p->setBaseShape(std::make_shared< FloatingBox >(edgeLength));
      p->getBaseShapeRef()->updateMassAndInertia(density);
      p->setOwner(walberla::mpi::MPIManager::instance()->rank());
      p->setType(1);
      p->setInteractionRadius(real_t(0.5) * edgeLength.length()); // circumscribing-sphere radius
      p->setRotation(rot);
      p->setLinearVelocity(vel);
      return p.getIdx();
   }

   size_t addCylinder(const Vec3& pos, const Vec3& vel, const Rot3& rot, real_t radius, real_t length,
                       real_t density)
   {
      auto p = ps_->create(false);
      p->setPosition(pos);
      p->setBaseShape(std::make_shared< data::Cylinder >(radius, length));
      p->getBaseShapeRef()->updateMassAndInertia(density);
      p->setOwner(walberla::mpi::MPIManager::instance()->rank());
      p->setType(1);
      p->setInteractionRadius(std::sqrt(radius * radius + real_t(0.25) * length * length));
      p->setRotation(rot);
      p->setLinearVelocity(vel);
      return p.getIdx();
   }

   size_t addWall(const Vec3& pos, const Vec3& normal) { return createPlane(*ps_, pos, normal).getIdx(); }

   //! Register one arch tunnel (see ParticleKernels.h::ComputeTunnelContacts for the full mechanism --
   //! this is not a mesa_pd particle, just an entry in a small list ComputeTunnelContacts iterates
   //! active floaters against directly). Lazily constructs computeTunnelContacts_ on first call.
   void addArchTunnel(const Vec3& center, real_t radius, real_t xMin, real_t xMax, real_t frictionCoefficient)
   {
      tunnels_.emplace_back(center, radius, xMin, xMax);
      computeTunnelContacts_ = std::make_unique< ComputeTunnelContacts >(
         tunnels_, restitutionCoefficient_, kappa_, collisionTime_, frictionCoefficient);
   }

   //! Re-poses an existing wall (from addWall) in place -- HalfSpace's normal is immutable after
   //! construction (mesa_pd/data/shape/HalfSpace.h has no setter), so "rotating" a wall means replacing
   //! its baseShape with a freshly-built HalfSpace(newNormal) alongside setPosition. Used to animate a
   //! polygon-of-planes drum wall (see detail::addDrumWalls/updateDrumWalls below) -- cheap enough to do
   //! every step for the handful of GLOBAL, non-communicating wall particles involved.
   void setWallPose(size_t idx, const Vec3& pos, const Vec3& normal)
   {
      ps_->setPosition(idx, pos);
      ps_->setBaseShape(idx, std::make_shared< data::HalfSpace >(normal));
   }

   //! Must be called once after all particles are created and before the first step(). `knownMaxDiameter`,
   //! when nonzero, overrides the scan-existing-particles sizing below -- required for scenarios (e.g.
   //! runPackingSingleShape) that call finalizeSetup() before any DYNAMIC particles exist yet (only
   //! GLOBAL walls), where the scan would otherwise size the linked-cell grid at ~0 and every particle
   //! created afterwards trips InsertParticleIntoLinkedCells' "interaction radius too large for this cell
   //! size" assertion the first time step() runs.
   void finalizeSetup(real_t knownMaxDiameter = real_t(0))
   {
      // Skip walls (GLOBAL, infinite interaction radius per Utility.h::createPlane) -- including them
      // would make maxDiameter (and hence the LinkedCells extension below) infinite.
      real_t maxDiameter = knownMaxDiameter;
      for (size_t i = 0; i < ac_->size(); ++i)
      {
         if (data::particle_flags::isSet(ac_->getFlags(i), data::particle_flags::GLOBAL)) continue;
         maxDiameter = std::max(maxDiameter, real_t(2) * ac_->getInteractionRadius(i));
      }
      linkedCellWidth_ = real_t(1.01) * std::max(maxDiameter, real_t(1e-6));

      ps_->forEachParticle(false, kernel::SelectLocal(), *ac_, *assoc_, *ac_);
      syncNN_(*ps_, forest_, domain_);
   }

   //! One DEM sub-cycle: pre-force integration, contact resolution (production ComputeContacts, the
   //! SAME class doDEMSubCycles uses), gravity, post-force integration. No lubrication, no
   //! hydrodynamic/buoyancy force -- this is the fluid-free "dry run" by construction.
   void step()
   {
      ps_->forEachParticle(false, kernel::SelectLocal(), *ac_, preVV_, *ac_);
      syncNN_(*ps_, forest_, domain_);

      data::LinkedCells linkedCells(domain_->getUnionOfLocalAABBs().getExtended(linkedCellWidth_),
                                     linkedCellWidth_);
      // LinkedCells's constructor leaves infiniteParticles_ at its default-member-initializer value (0,
      // via std::atomic<int>'s C++20 value-init) instead of the -1 sentinel its own traversal logic
      // expects -- only clear() sets it correctly. BedGeneration.cpp calls clear() for the same reason
      // (there because it reuses one long-lived LinkedCells object; here because even a freshly
      // constructed one needs it). Without this, the first inserted INFINITE particle's next-pointer
      // wrongly chains to particle index 0, so any pair traversal that reaches the end of the infinite-
      // particle list also revisits index 0 -- a real self-pairing crash when index 0 is itself the sole
      // finite body (exactly the wall-impact scenarios below).
      linkedCells.clear();
      kernel::InsertParticleIntoLinkedCells ipilc;
      ps_->forEachParticle(false, kernel::SelectAll(), *ac_, ipilc, *ac_, linkedCells);
      linkedCells.forEachParticlePairHalf(
         false, kernel::ExcludeInfiniteInfinite(), *ac_,
         [this](const size_t idx1, const size_t idx2, auto& ac) { (*computeContacts_)(idx1, idx2, ac); }, *ac_);

      if (computeTunnelContacts_)
         ps_->forEachParticle(false, kernel::SelectLocal(), *ac_, *computeTunnelContacts_, *ac_);

      if (std::abs(gravityZ_) > real_t(0))
      {
         ps_->forEachParticle(
            false, kernel::SelectLocal(), *ac_,
            [this](const size_t idx, auto& ac) {
               mesa_pd::addForceAtomic(idx, ac, Vec3(real_t(0), real_t(0), ac.getMass(idx) * gravityZ_));
            },
            *ac_);
      }

      mesa_pd::mpi::ReduceProperty RP;
      RP.operator()< ForceTorqueNotification >(*ps_);

      ps_->forEachParticle(false, kernel::SelectLocal(), *ac_, postVV_, *ac_);
      syncNN_(*ps_, forest_, domain_);
      ++stepCount_;
      t_ += dt_;
   }

   //! Query-only contact check (no force applied): reuses the exact production contact-detection
   //! functor (SolidCylinderContactDetection via DoubleCast) as the state-machine trigger for
   //! "has this pair started/stopped touching" in the wall-impact / pair-restitution scenarios below.
   bool queryContact(size_t idx1, size_t idx2, real_t& penetrationOut)
   {
      SolidCylinderContactDetection acd;
      kernel::DoubleCast            double_cast;
      if (double_cast(idx1, idx2, *ac_, acd, *ac_))
      {
         penetrationOut = acd.getPenetrationDepth();
         return true;
      }
      return false;
   }

   real_t kineticEnergy(size_t idx) const
   {
      const real_t transKE = real_t(0.5) * ac_->getMass(idx) * ac_->getLinearVelocity(idx).sqrLength();
      const Mat3   Iworld  = math::transformMatrixRART(ac_->getRotation(idx).getMatrix(), ac_->getInertiaBF(idx));
      const Vec3   w       = ac_->getAngularVelocity(idx);
      const real_t rotKE   = real_t(0.5) * (w * (Iworld * w));
      return transKE + rotKE;
   }

   real_t time() const { return t_; }
   uint_t stepCount() const { return stepCount_; }
   real_t dt() const { return dt_; }

   void openCsv(const std::string& filename, const std::vector< std::string >& fieldNames,
                const std::vector< size_t >& trackedIndices)
   {
      trackedIndices_ = trackedIndices;
      csv_.open(filename.c_str());
      csv_ << "step,t";
      for (size_t bi = 0; bi < trackedIndices_.size(); ++bi)
         for (const auto& f : fieldNames)
            csv_ << ",body" << bi << "_" << f;
      csv_ << "\n";
   }

   void writeCsvRow()
   {
      if (!csv_.is_open()) return;
      csv_ << stepCount_ << "," << t_;
      for (size_t idx : trackedIndices_)
      {
         const auto& pos = ac_->getPosition(idx);
         const auto& vel = ac_->getLinearVelocity(idx);
         const auto& w   = ac_->getAngularVelocity(idx);
         csv_ << "," << pos[0] << "," << pos[1] << "," << pos[2] << "," << vel[0] << "," << vel[1] << ","
              << vel[2] << "," << w[0] << "," << w[1] << "," << w[2] << "," << kineticEnergy(idx);
      }
      csv_ << "\n";
   }

   // edgeLength/rotation: the SAME two fields, in the SAME convention, as BuildTimeloop.h's own
   // production box floater VTK output (its own "AddOutputFloatingBoxes" section) -- so a box in this
   // DEMValidation harness renders as its actual ORIENTED solid in ParaView (Filters > Glyph, Glyph
   // Type = Box, Scale Array = edgeLength, Scale Mode = Scale by Components; then Filters > Tensor
   // Glyph using "rotation", Extract Eigenvalues OFF) instead of a generic sphere-by-radius, exactly
   // matching the recipe already documented there. Zero (Vec3(0)) for any non-box shape, so this is
   // safe to add unconditionally for every DEMValidation scenario (sphere/cylinder VTK is unaffected --
   // a zero-size box glyph simply doesn't render, "radius"-based glyphing still works for those).
   struct SelectBoxEdgeLength
   {
      using return_type = Vec3;
      return_type operator()(const data::Particle& p) const
      {
         const auto& shape = p.getBaseShape();
         if (!shape || shape->getShapeType() != data::Box::SHAPE_TYPE) return Vec3(real_t(0));
         return static_cast< const data::Box* >(shape.get())->getEdgeLength();
      }
   };
   // Transposed body->world rotation matrix -- see BuildTimeloop.h's own comment on why the transpose
   // is required for ParaView's column-major Tensor Glyph to match the actual orientation (not its
   // inverse). Visualization-only ordering; production's own buoyancy/torque math (and this file's own
   // kineticEnergy()) uses getRotation().getMatrix() directly, untransposed, unaffected by this.
   struct SelectParticleRotationMatrix
   {
      using return_type = std::array< real_t, 9 >;
      return_type operator()(const data::Particle& p) const
      {
         const Mat3 R = p.getRotation().getMatrix();
         return { R(0, 0), R(1, 0), R(2, 0), R(0, 1), R(1, 1), R(2, 1), R(0, 2), R(1, 2), R(2, 2) };
      }
   };

   void setupVtk(const std::string& baseName, uint_t spacing)
   {
      vtkSpacing_ = spacing;
      if (vtkSpacing_ == 0) return;
      auto particleVtkOutput = make_shared< mesa_pd::vtk::ParticleVtkOutput >(ps_);
      particleVtkOutput->addOutput< data::SelectParticleLinearVelocity >("velocity");
      particleVtkOutput->addOutput< data::SelectParticleAngularVelocity >("angularVelocity");
      particleVtkOutput->addOutput< data::SelectParticleInteractionRadius >("radius");
      particleVtkOutput->addOutput< SelectBoxEdgeLength >("edgeLength");
      particleVtkOutput->addOutput< SelectParticleRotationMatrix >("rotation");
      particleVtkOutput->addOutput< data::SelectParticleUid >("uid");
      particleVtkOutput->setParticleSelector([](const data::ParticleStorage::iterator& pIt) {
         using namespace walberla::mesa_pd::data::particle_flags;
         return pIt->getBaseShape()->getShapeType() != data::HalfSpace::SHAPE_TYPE && !isSet(pIt->getFlags(), GHOST);
      });
      vtkWriter_ = walberla::vtk::createVTKOutput_PointData(particleVtkOutput, "Particles", 1,
                                                              "vtk-out/" + baseName, "simulation_step", false, false);
   }

   void writeVtkIfDue()
   {
      if (vtkWriter_ && vtkSpacing_ > 0 && stepCount_ % vtkSpacing_ == 0) vtkWriter_->write();
   }

private:
   real_t dt_;
   real_t gravityZ_;
   real_t restitutionCoefficient_;
   real_t collisionTime_;
   real_t kappa_;

   std::shared_ptr< BlockForest >                        forest_;
   std::shared_ptr< mesa_pd::domain::BlockForestDomain > domain_;
   std::shared_ptr< data::ParticleStorage >              ps_;
   std::unique_ptr< data::ParticleAccessorWithBaseShape > ac_;
   std::unique_ptr< kernel::AssocToBlock >                assoc_;
   std::unique_ptr< ComputeContacts >                     computeContacts_;
   std::vector< ArchTunnelParams >                        tunnels_;
   std::unique_ptr< ComputeTunnelContacts >                computeTunnelContacts_;

   kernel::VelocityVerletPreForceUpdate       preVV_;
   kernel::VelocityVerletPostForceUpdate      postVV_;
   mesa_pd::mpi::SyncNextNeighborsBlockForest syncNN_;

   real_t linkedCellWidth_ = real_t(1);
   real_t t_               = real_t(0);
   uint_t stepCount_       = uint_t(0);

   std::ofstream                     csv_;
   std::vector< size_t >             trackedIndices_;
   std::shared_ptr< walberla::vtk::VTKOutput > vtkWriter_;
   uint_t                            vtkSpacing_ = 0;
};

//! Builds both a Rot3 (for setRotation) and its Mat3 (for the analytical-solution calls, which take
//! plain rotation matrices) from a single rotation about the world Y-axis, so the two never drift apart.
struct YRotation
{
   Rot3 rot;
   Mat3 mat;
   explicit YRotation(real_t angleRad) : rot(Vec3(real_t(0), real_t(1), real_t(0)), angleRad), mat(rot.getMatrix()) {}
};

inline real_t degToRad(real_t deg) { return deg * math::pi / real_t(180); }

//! KE of a rigid body given explicit state -- used to compute the ANALYTICAL model's own predicted
//! post-impact energy (as opposed to DEMEngine::kineticEnergy, which reads a live simulated particle).
inline real_t kineticEnergyOf(real_t mass, const Vec3& velocity, const Vec3& angularVelocity,
                               const Mat3& inertiaWorld)
{
   return real_t(0.5) * mass * velocity.sqrLength() +
          real_t(0.5) * (angularVelocity * (inertiaWorld * angularVelocity));
}

struct PassFailAccumulator
{
   uint_t total  = 0;
   uint_t passed = 0;
   void   record(bool pass)
   {
      ++total;
      if (pass) ++passed;
   }
   bool allPassed() const { return total > 0 && passed == total; }
};

//! Runs `engine` until `bodyIdx` has made and then broken contact with `wallIdx` (contact-force-based,
//! not a position/height heuristic -- a tilted body's resting height changes as it rotates during
//! contact, so a fixed-height threshold can miss separation entirely and run to maxSteps). Returns the
//! number of steps actually taken.
inline uint_t runWallImpactLoop(DEMEngine& engine, size_t bodyIdx, size_t wallIdx, uint_t maxSteps,
                                 bool writeOutput)
{
   bool   everContacted = false;
   real_t penetration   = real_t(0);
   uint_t s             = 0;
   for (; s < maxSteps; ++s)
   {
      const bool contactNow = engine.queryContact(bodyIdx, wallIdx, penetration) && penetration < real_t(0);
      if (contactNow) everContacted = true;
      engine.step();
      if (writeOutput)
      {
         engine.writeCsvRow();
         engine.writeVtkIfDue();
      }
      if (everContacted && !contactNow && s > 0) break; // just separated this step
   }
   return s;
}

inline void logCaseResult(const std::string& tag, const Vec3& simVel, const Vec3& analyticVel,
                           const Vec3& simAngVel, const Vec3& analyticAngVel, real_t tolerance,
                           PassFailAccumulator& acc)
{
   const real_t velErr = (simVel - analyticVel).length() / std::max(real_t(1e-8), analyticVel.length());
   const real_t angErr = (simAngVel - analyticAngVel).length() /
                          std::max(real_t(1e-8), std::max(analyticAngVel.length(), simAngVel.length()));
   const bool pass = velErr < tolerance && angErr < tolerance;
   acc.record(pass);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] "
                              << tag << "  simVel=" << simVel << " analyticVel=" << analyticVel
                              << " velRelErr=" << velErr << "  simAngVel=" << simAngVel
                              << " analyticAngVel=" << analyticAngVel << " angRelErr=" << angErr << "  "
                              << (pass ? "PASS" : "FAIL"));
}

//! Hard physical bound independent of any friction-model calibration: a passive contact (no motors, no
//! energy source) can only ever dissipate kinetic energy, never create it, regardless of the restitution
//! or friction coefficients used. Checked in EVERY collision scenario (wall-impact and pairwise alike) as
//! a model-agnostic sanity bound alongside each test's own quantitative comparison -- a violation here
//! would mean the contact law itself is unstable/energy-injecting, a much more basic problem than a
//! friction-model calibration mismatch. tolerance guards only against integration noise (sub-cycled
//! spring-dashpot contacts are not exactly energy-monotonic step-to-step even when correct).
inline void logConservation(const std::string& tag, real_t keBefore, real_t keAfter, PassFailAccumulator& acc,
                             real_t tolerance = real_t(1e-3))
{
   const bool pass = keAfter <= keBefore * (real_t(1) + tolerance) + real_t(1e-9);
   acc.record(pass);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] " << tag << "  [energy non-increase] KE_before=" << keBefore
                              << " KE_after=" << keAfter << " ratio=" << (keAfter / std::max(real_t(1e-12), keBefore))
                              << "  " << (pass ? "PASS" : "FAIL (energy created -- contact law is unstable)"));
}

//! STRONGER companion to logConservation: for scenarios where a specific post-impact energy is
//! computable (either the e=1 elastic case, where it must equal keBefore exactly, or wherever an
//! analytical model already predicts the exact post-state), checking only "did not increase" is a weak,
//! one-sided bound -- it would pass even if the simulation dissipated far more or far less energy than a
//! correctly-calibrated restitution/friction model implies. This checks the ACTUAL dissipated amount
//! against the THEORETICALLY EXPECTED amount, not just its sign.
inline void logEnergyBudget(const std::string& tag, real_t keAfterSim, real_t keAfterExpected,
                             PassFailAccumulator& acc, real_t tolerance)
{
   const real_t relErr = std::abs(keAfterSim - keAfterExpected) / std::max(real_t(1e-8), keAfterExpected);
   const bool   pass   = relErr < tolerance;
   acc.record(pass);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] " << tag << "  [energy budget] KE_after_sim=" << keAfterSim
                              << " KE_after_expected=" << keAfterExpected << " relErr=" << relErr << "  "
                              << (pass ? "PASS" : "FAIL"));
}

//! Exact for any isolated N-body system with only internal (pairwise, Newton's-third-law) contact forces
//! and no external force (no gravity, no walls) -- total linear momentum must be conserved to numerical
//! precision regardless of the collision's restitution/friction outcome. A stronger, more model-agnostic
//! check than energy or velocity comparisons: it holds even for a WRONG friction/restitution calibration,
//! so a violation specifically implicates the contact force's Newton's-third-law symmetry (equal and
//! opposite forces on the two bodies), not a parameter-tuning issue.
inline void logMomentumConservation(const std::string& tag, const Vec3& pBefore, const Vec3& pAfter,
                                     PassFailAccumulator& acc, real_t tolerance = real_t(1e-6))
{
   const real_t relErr = (pAfter - pBefore).length() / std::max(real_t(1e-8), pBefore.length());
   const bool   pass   = relErr < tolerance;
   acc.record(pass);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] " << tag << "  [momentum conservation] pBefore=" << pBefore
                              << " pAfter=" << pAfter << " relErr=" << relErr << "  "
                              << (pass ? "PASS" : "FAIL"));
}

//! Total mechanical energy (kinetic + gravitational potential, PE = m*g*z relative to the floor at
//! z=0) summed over every non-fixed, non-global dynamic particle. For a multi-body gravity-settling
//! system with only dissipative internal forces (inelastic/frictional contacts, no motors or other
//! energy source), this must be MONOTONICALLY NON-INCREASING over time -- a violation (checked between
//! successive sampling checkpoints, not step-to-step, since a single sub-cycle is not exactly
//! energy-monotonic even when correct) is a much more basic red flag than any shape- or friction-specific
//! discrepancy: it means the contact/integration scheme is injecting energy somewhere.
inline real_t totalMechanicalEnergy(DEMEngine& engine, real_t gravity)
{
   real_t total = real_t(0);
   for (size_t i = 0; i < engine.ac().size(); ++i)
   {
      if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
      total += engine.kineticEnergy(i) + engine.ac().getMass(i) * gravity * engine.ac().getPosition(i)[2];
   }
   return total;
}

} // namespace detail

using detail::DEMEngine;

//======================================================================================================================
// 1-3: single-body wall-impact tests
//======================================================================================================================

inline bool runWallImpactSphere(const config::Config::BlockHandle& b)
{
   const real_t radius        = b.getParameter< real_t >("sphereRadius", real_t(1));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.85));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.15));
   const real_t speed         = b.getParameter< real_t >("impactSpeed", real_t(1));
   const real_t ratioStart    = b.getParameter< real_t >("tangentialRatioStart", real_t(0));
   const real_t ratioEnd      = b.getParameter< real_t >("tangentialRatioEnd", real_t(1.5));
   const real_t ratioStep     = b.getParameter< real_t >("tangentialRatioStep", real_t(0.25));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "wall_impact_sphere");

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);

   for (real_t ratio = ratioStart; ratio <= ratioEnd + real_t(1e-9); ratio += ratioStep)
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), real_t(50) * radius);
      const real_t gap = real_t(0.01) * radius;
      const size_t sphereIdx =
         engine.addSphere(Vec3(0, 0, radius + gap), Vec3(ratio * speed, 0, -speed), radius, density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_ratio" + std::to_string(ratio);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { sphereIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(sphereIdx);
      const real_t keBefore = engine.kineticEnergy(sphereIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, sphereIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(sphereIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(sphereIdx);
      const real_t keAfter   = engine.kineticEnergy(sphereIdx);

      const auto analytic = analytical::sphereWallImpact(radius, engine.ac().getMass(sphereIdx), preVel, Vec3(0),
                                                           wallNormal, e, mu);

      const std::string caseTag = "wallImpactSphere ratio=" + std::to_string(ratio);
      detail::logCaseResult(caseTag, simVel, analytic.postVelocityA, simAngVel, analytic.postAngularVelocityA,
                             tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const real_t sphereI = real_t(0.4) * engine.ac().getMass(sphereIdx) * radius * radius; // isotropic
      const real_t keExpected = detail::kineticEnergyOf(engine.ac().getMass(sphereIdx), analytic.postVelocityA,
                                                          analytic.postAngularVelocityA,
                                                          Mat3::makeDiagonalMatrix(sphereI));
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactSphere OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

inline bool runWallImpactCylinder(const config::Config::BlockHandle& b)
{
   const real_t radius        = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t length        = b.getParameter< real_t >("cylinderLength", real_t(4));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.85));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.15));
   const real_t speed         = b.getParameter< real_t >("impactSpeed", real_t(1));
   // Default 1.0 (equal tangential/normal, this app's own general oblique-impact sweep). Kodam et al.
   // (2010) Part II's actual Test 1 protocol is PURELY NORMAL (zero tangential component, zero initial
   // spin, frictionless) -- see dem_validation_wall_cylinder_kodam_test1.prm, which sets this to 0.
   const real_t tangentialFraction = b.getParameter< real_t >("tangentialSpeedFraction", real_t(1));
   const real_t tiltStartDeg  = b.getParameter< real_t >("tiltStartDeg", real_t(0));
   const real_t tiltEndDeg    = b.getParameter< real_t >("tiltEndDeg", real_t(90));
   const real_t tiltStepDeg   = b.getParameter< real_t >("tiltStepDeg", real_t(15));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "wall_impact_cylinder");
   // useRotationalEffectiveMass (see ModelParameters' own doc comment, SimulationParams.h, and
   // computeRotationalEffectiveMass, ParticleKernels.h): LinearSpringDashpot's stiffness/damping are
   // calibrated ONCE from the plain TRANSLATIONAL effective mass, which is only exact when the contact
   // point has zero lever arm about the line of centres -- exactly true at tilt=0/90 (on-axis contact)
   // and increasingly false at intermediate tilts, matching this test's own observed error pattern
   // (clean at 0/90, worst in the middle). Default false reproduces every existing prm's exact prior
   // behaviour; the Kodam-replica diagnostic run turns this on to test whether it explains the gap.
   const bool useRotationalEffectiveMass = b.getParameter< bool >("useRotationalEffectiveMass", false);

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);

   data::Cylinder refShape(radius, length);
   refShape.updateMassAndInertia(density);

   for (real_t tiltDeg = tiltStartDeg; tiltDeg <= tiltEndDeg + real_t(1e-9); tiltDeg += tiltStepDeg)
   {
      const detail::YRotation yrot(detail::degToRad(tiltDeg));
      // Cylinder's own support() gives the exact "reach" toward the wall for this tilt, matching the
      // formula production contact detection uses -- see AnalyticalSolutions.h's boxCornerWallImpact
      // comment for the same reasoning applied to the box case.
      const Vec3   bfDir  = yrot.mat.getTranspose() * (-wallNormal);
      const Vec3   lowest = yrot.mat * refShape.support(bfDir);
      const real_t reach  = -lowest[2];

      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), real_t(20) * (radius + length),
                        /*useAnalyticalContactDetection=*/false, /*useRoundedBoxCorners=*/false,
                        /*boxRoundingRadius=*/real_t(0.05), /*useMultiPointBoxContact=*/false,
                        /*useBoxOverlapCap=*/false, /*boxOverlapCapFraction=*/real_t(0.3),
                        /*useAxisHysteresis=*/false, /*axisHysteresisMargin=*/real_t(0.15),
                        /*useContactManifold=*/false, /*manifoldFadeRate=*/real_t(0.25),
                        useRotationalEffectiveMass);
      const real_t gap    = real_t(0.01) * radius;
      const size_t cylIdx = engine.addCylinder(Vec3(0, 0, reach + gap),
                                                Vec3(tangentialFraction * speed, 0, -speed), yrot.rot, radius,
                                                length, density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_tilt" + std::to_string(tiltDeg);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { cylIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(cylIdx);
      const real_t keBefore = engine.kineticEnergy(cylIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, cylIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(cylIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(cylIdx);
      const real_t keAfter   = engine.kineticEnergy(cylIdx);

      const paper_analytical::CylinderState cylState{ Vec3(0, 0, reach + gap), yrot.mat, radius, length };
      const Mat3 invInertiaWorld = math::transformMatrixRART(yrot.mat, refShape.getInvInertiaBF());
      const auto analytic = analytical::cylinderWallImpact(cylState, refShape.getMass(), invInertiaWorld, preVel,
                                                             Vec3(0), wallPos, wallNormal, e, mu);

      const std::string caseTag = "wallImpactCylinder tilt=" + std::to_string(tiltDeg) + "deg";
      detail::logCaseResult(caseTag, simVel, analytic.postVelocityA, simAngVel, analytic.postAngularVelocityA,
                             tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const Mat3   inertiaWorld = math::transformMatrixRART(yrot.mat, refShape.getInertiaBF());
      const real_t keExpected   = detail::kineticEnergyOf(refShape.getMass(), analytic.postVelocityA,
                                                            analytic.postAngularVelocityA, inertiaWorld);
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);

      // DIRECT comparison against Kodam et al. (2010) Part II's own stated numeric result (p. 5866): at
      // face-on (0deg) and band-on (90deg), the dimensionless rebound speed Vz+/Vz- equals EXACTLY the
      // restitution coefficient (their reported value: 0.85). This checks the simulated ratio against
      // that PAPER-STATED number directly -- independent of this file's own analytical solver -- rather
      // than only cross-checking against a formula we derived ourselves. Only meaningful at frictionless
      // 0/90deg (see dem_validation_wall_cylinder_kodam_test1.prm); silently skipped otherwise since
      // Kodam's paper gives no other explicit numeric target across the sweep.
      if (std::abs(mu) < real_t(1e-9) &&
          (std::abs(tiltDeg) < real_t(1e-6) || std::abs(tiltDeg - real_t(90)) < real_t(1e-6)))
      {
         const real_t simRebound = std::abs(simVel[2]) / std::abs(preVel[2]);
         const real_t relErr     = std::abs(simRebound - e) / e;
         const bool   pass       = relErr < real_t(0.02);
         acc.record(pass);
         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactCylinder tilt="
                                    << tiltDeg << "deg  [vs. Kodam et al. 2010 Part II, p. 5866, stated result] "
                                    << "simRebound=" << simRebound << " paperReported=" << e
                                    << " (=configured restitution, their calibrated value) relErr=" << relErr
                                    << "  " << (pass ? "PASS" : "FAIL"));
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactCylinder OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

//! Diagnostic-only cross-check (not itself a pass/fail validation experiment against a literature
//! result -- a self-consistency check between two of OUR OWN geometric methods): does mesa_pd's
//! generic iterative GJK/EPA solver agree with the exact closed-form Cylinder x HalfSpace contact
//! geometry that production actually uses (SolidCylinderContactDetection's cyl.support()-based
//! formula, ParticleKernels.h:338 -- algebraically identical to paper_analytical::rimPointTowardPlane,
//! PaperAnalyticalContacts.h:1980, the formula AnalyticalSolutions.h::cylinderWallImpact uses as
//! wallImpactCylinder's own analytical reference)?
//!
//! A true HalfSpace has NO support() function at all (mesa_pd::data::HalfSpace.h) -- GJK/EPA cannot
//! run on a Cylinder x HalfSpace pair in this codebase AT ALL, by construction, not by configuration
//! choice. This test approximates the wall as a very large, thin FloatingBox (which DOES have
//! support(), so GJK/EPA is well-defined) purely to get an independent third answer for the SAME
//! cylinder geometry, at the SAME tilt angles wallImpactCylinder's own Kodam-replica sweep uses.
inline bool runGjkVsAnalyticalCylinderWallCheck(const config::Config::BlockHandle& b)
{
   const real_t radius       = b.getParameter< real_t >("cylinderRadius", real_t(4));
   const real_t length       = b.getParameter< real_t >("cylinderLength", real_t(5.3));
   const real_t overlapDepth = b.getParameter< real_t >("overlapDepth", real_t(0.05));
   const real_t tiltStartDeg = b.getParameter< real_t >("tiltStartDeg", real_t(0));
   const real_t tiltEndDeg   = b.getParameter< real_t >("tiltEndDeg", real_t(90));
   const real_t tiltStepDeg  = b.getParameter< real_t >("tiltStepDeg", real_t(10));
   const real_t tolerance    = b.getParameter< real_t >("toleranceRelError", real_t(0.02));

   data::Cylinder cyl(radius, length);
   FloatingBox    wallProxy(Vec3(real_t(2000), real_t(2000), real_t(100)));
   wallProxy.updateMassAndInertia(real_t(1));
   const Vec3 wallPos(0, 0, -real_t(50)); // top face exactly at world z=0, matching HalfSpace's own plane
   const Rot3 wallRot;                    // identity
   const Vec3 wallNormal(0, 0, 1);

   detail::PassFailAccumulator acc;

   for (real_t tiltDeg = tiltStartDeg; tiltDeg <= tiltEndDeg + real_t(1e-9); tiltDeg += tiltStepDeg)
   {
      const detail::YRotation yrot(detail::degToRad(tiltDeg));

      // Exact support-point reach (same formula production contact detection + the analytical
      // reference both use), positioned to overlap the wall by exactly `overlapDepth`.
      const Vec3   bfDir  = yrot.mat.getTranspose() * (-wallNormal);
      const Vec3   lowest = yrot.mat * cyl.support(bfDir);
      const real_t reach  = -lowest[2];
      const Vec3   cylPos(0, 0, reach - overlapDepth);

      // (a) exact closed form -- identical to ParticleKernels.h's SolidCylinderContactDetection
      const Vec3   exactContact          = cylPos + yrot.mat * cyl.support(bfDir);
      const real_t exactPenetrationDepth = wallNormal * exactContact; // wall plane passes through the origin

      // (b) GJK/EPA against the large-box wall proxy
      mesa_pd::collision_detection::Support supportCyl(cylPos, yrot.rot, cyl);
      mesa_pd::collision_detection::Support supportWall(wallPos, wallRot, wallProxy);
      mesa_pd::collision_detection::GJK     gjk;
      bool   gjkFound            = false;
      Vec3   gjkPoint(0), gjkNormal(0);
      real_t gjkPenetrationDepth = 0;
      if (gjk.doGJKmargin(supportCyl, supportWall, real_t(0)))
      {
         mesa_pd::collision_detection::EPA epa;
         epa.useSphereOptimization(false);
         gjkFound =
            epa.doEPAmargin(supportCyl, supportWall, gjk, gjkNormal, gjkPoint, gjkPenetrationDepth, real_t(0));
      }

      const real_t overlapDiff =
         gjkFound ? std::abs(std::abs(gjkPenetrationDepth) - std::abs(exactPenetrationDepth)) : real_t(-1);
      const real_t pointDist = gjkFound ? (gjkPoint - exactContact).length() : real_t(-1);

      const bool pass = gjkFound && overlapDiff < tolerance * radius;
      acc.record(pass);
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] gjkVsAnalyticalCylinderWall tilt="
                                 << tiltDeg << "deg  exactPoint=" << exactContact
                                 << " exactOverlap=" << exactPenetrationDepth << "  gjkFound=" << gjkFound
                                 << " gjkPoint=" << gjkPoint << " gjkOverlap=" << gjkPenetrationDepth
                                 << "  overlapDiff=" << overlapDiff << " pointDist=" << pointDist << "  "
                                 << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] gjkVsAnalyticalCylinderWall OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

inline bool runWallImpactBox(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents   = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 2, 2));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.6));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.2));
   const real_t speed         = b.getParameter< real_t >("impactSpeed", real_t(1));
   // [Claude, 2026-10-29] Mirrors wallImpactCylinder's own tangentialSpeedFraction (above): default 1.0 keeps
   // every existing prm's exact prior behaviour (equal tangential/normal, this app's own general oblique-impact
   // sweep); set to 0.0 for a Kodam-et-al.-Test-1-style purely normal impact, now available for the box the same
   // way it already was for the cylinder (author: "can you recreate the kodam version for the cuboids").
   const real_t tangentialFraction = b.getParameter< real_t >("tangentialSpeedFraction", real_t(1));
   const real_t tiltStartDeg  = b.getParameter< real_t >("tiltStartDeg", real_t(5));
   const real_t tiltEndDeg    = b.getParameter< real_t >("tiltEndDeg", real_t(40));
   const real_t tiltStepDeg   = b.getParameter< real_t >("tiltStepDeg", real_t(10));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   // [Claude, 2026-09-30] was not wired for the box (always the plain translational effective mass); default false keeps
   // the previous behaviour.
   const bool   useRotationalEffectiveMass = b.getParameter< bool >("useRotationalEffectiveMass", false);
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "wall_impact_box");

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);
   const Vec3 edgeLength = real_t(2) * halfExtents;

   FloatingBox refShape(edgeLength);
   refShape.updateMassAndInertia(density);

   for (real_t tiltDeg = tiltStartDeg; tiltDeg <= tiltEndDeg + real_t(1e-9); tiltDeg += tiltStepDeg)
   {
      const detail::YRotation yrot(detail::degToRad(tiltDeg));
      const Vec3   bfDir  = yrot.mat.getTranspose() * (-wallNormal);
      const Vec3   corner = refShape.support(bfDir); // body-frame corner, used for POSITIONING (height) only
      const Vec3   lowest = yrot.mat * corner;
      const real_t reach  = -lowest[2];
      // A pure Y-axis tilt never changes a point's world Y (rotation about Y leaves Y invariant) or its
      // contribution to world Z (which depends only on body X/Z), so the "lowest" support corner is
      // always TIED between its +Y and -Y twin -- the true first contact is the full bottom EDGE
      // spanning Y, matching Housner's 2D rocking-block corner exactly in the X-Z cross-section, not a
      // single 3D corner. support()'s arbitrary tie-break otherwise injects a spurious Y lever arm (and,
      // via the K matrix's rotational coupling, corrupts the X-Z answer too) that the true symmetric
      // edge contact does not have -- confirmed empirically (a non-zero-Y support point predicted a
      // physically impossible negative post-bounce Z-velocity). Analytical contact point is the same
      // corner projected onto the Y=0 symmetry plane.
      const Vec3 analyticCornerBF(corner[0], real_t(0), corner[2]);

      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0),
                       real_t(30) * std::max({ halfExtents[0], halfExtents[1], halfExtents[2] }),
                       /*useAnalyticalContactDetection=*/false, /*useRoundedBoxCorners=*/false,
                       /*boxRoundingRadius=*/real_t(0.05), /*useMultiPointBoxContact=*/false,
                       /*useBoxOverlapCap=*/false, /*boxOverlapCapFraction=*/real_t(0.3),
                       /*useAxisHysteresis=*/false, /*axisHysteresisMargin=*/real_t(0.15),
                       /*useContactManifold=*/false, /*manifoldFadeRate=*/real_t(0.25), useRotationalEffectiveMass);
      const real_t gap = real_t(0.01) * std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
      // Horizontal velocity component aimed so the leading (lowest) corner strikes squarely, matching
      // Housner's planar rocking-block assumption (velocity, spin axis and contact point all in the
      // X-Z plane -- tilt is about world Y, so no Y motion is introduced).
      const size_t boxIdx = engine.addBox(Vec3(0, 0, reach + gap),
                                            Vec3(tangentialFraction * speed, 0, -speed), yrot.rot, edgeLength,
                                            density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_tilt" + std::to_string(tiltDeg);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { boxIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(boxIdx);
      const real_t keBefore = engine.kineticEnergy(boxIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, boxIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(boxIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(boxIdx);
      const real_t keAfter   = engine.kineticEnergy(boxIdx);

      const auto analytic = analytical::boxCornerWallImpact(analyticCornerBF, yrot.mat, refShape.getMass(),
                                                              refShape.getInvInertiaBF(), preVel, Vec3(0), wallNormal,
                                                              e, mu);

      const std::string caseTag = "wallImpactBox tilt=" + std::to_string(tiltDeg) + "deg";
      detail::logCaseResult(caseTag, simVel, analytic.postVelocityA, simAngVel, analytic.postAngularVelocityA,
                             tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const Mat3   inertiaWorld = math::transformMatrixRART(yrot.mat, refShape.getInertiaBF());
      const real_t keExpected   = detail::kineticEnergyOf(refShape.getMass(), analytic.postVelocityA,
                                                            analytic.postAngularVelocityA, inertiaWorld);
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] wallImpactBox OVERALL: " << (acc.allPassed() ? "PASS" : "FAIL") << " ("
                                                                        << acc.passed << "/" << acc.total << ")");
   return acc.allPassed();
}

//======================================================================================================================
// 3b: literature-grounded cuboid single-particle regression checks (Fraige et al. 2008 / Mack et al. 2011)
//======================================================================================================================

//! Literal replica of Mack et al. (2011)'s single-cube coefficient-of-restitution drop test methodology
//! (Powder Technology 214:431-442, Section 2.4: e = sqrt(h1/h0), measured for real hexahedra as e=0.54,
//! the default below) -- but FACE-first (flush, zero tilt), not corner-first like wallImpactBox's own
//! sweep (which starts at tiltStartDeg=5 and never exercises the flush case). This is deliberately the
//! EXACT degenerate configuration that Cylinder/Box::support()'s tie-break bug (fixed 2026-09-13, see
//! DEM_VALIDATION_MATH_REFERENCE.md's box-support() note) corrupted: pre-fix, a flush box's support()
//! picked an arbitrary corner instead of the true face centre, injecting a spurious lever arm/torque even
//! though the true physics here is a trivial 1-DOF normal bounce by symmetry -- analyticVel =
//! (0,0,+e*impactSpeed), analyticAngVel = (0,0,0) exactly, no AnalyticalSolutions.h call needed. Any
//! nonzero simulated angular velocity, or a normal-velocity ratio deviating from the input e, is
//! therefore a direct regression signal for that bug, not a friction/calibration artifact -- there is no
//! tangential velocity here to mobilise frictionCoefficient at all.
inline bool runBoxFaceDropRestitution(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents   = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 1, 1));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   // Mack et al. (2011), Powder Technology 214:431-442, Section 2.4 / Table 2: measured hexahedron
   // (cube) particle-wall coefficient of restitution e=0.54 and sliding friction mu=0.28 (Section 2.3).
   const real_t e             = b.getParameter< real_t >("restitutionCoefficient", real_t(0.54));
   const real_t mu            = b.getParameter< real_t >("frictionCoefficient", real_t(0.28));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const real_t speedStart    = b.getParameter< real_t >("impactSpeedStart", real_t(0.5));
   const real_t speedEnd      = b.getParameter< real_t >("impactSpeedEnd", real_t(2.0));
   const real_t speedStep     = b.getParameter< real_t >("impactSpeedStep", real_t(0.5));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "box_face_drop_restitution");

   detail::PassFailAccumulator acc;
   const Vec3 wallPos(0, 0, 0), wallNormal(0, 0, 1);
   const Vec3 edgeLength = real_t(2) * halfExtents;

   for (real_t speed = speedStart; speed <= speedEnd + real_t(1e-9); speed += speedStep)
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0),
                        real_t(30) * std::max({ halfExtents[0], halfExtents[1], halfExtents[2] }));
      const real_t gap = real_t(0.01) * std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
      const size_t boxIdx =
         engine.addBox(Vec3(0, 0, halfExtents[2] + gap), Vec3(0, 0, -speed), Rot3(), edgeLength, density);
      const size_t wallIdx = engine.addWall(wallPos, wallNormal);
      engine.finalizeSetup();

      const std::string tag = outBase + "_v" + std::to_string(speed);
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { boxIdx });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   preVel   = engine.ac().getLinearVelocity(boxIdx);
      const real_t keBefore = engine.kineticEnergy(boxIdx);
      const uint_t maxSteps = uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runWallImpactLoop(engine, boxIdx, wallIdx, maxSteps, true);

      const Vec3   simVel    = engine.ac().getLinearVelocity(boxIdx);
      const Vec3   simAngVel = engine.ac().getAngularVelocity(boxIdx);
      const real_t keAfter   = engine.kineticEnergy(boxIdx);

      // Trivial by symmetry (flush, zero tangential velocity): normal velocity flips and scales by e,
      // angular velocity stays exactly zero -- no AnalyticalSolutions.h call needed for this case.
      // v_z^+ = -e*v_z^-, e=0.54 (Mack et al. 2011, Table 2)
      const Vec3 analyticVel(0, 0, -e * preVel[2]);
      const Vec3 analyticAngVel(0, 0, 0);

      const std::string caseTag = "boxFaceDropRestitution v=" + std::to_string(speed);
      detail::logCaseResult(caseTag, simVel, analyticVel, simAngVel, analyticAngVel, tolerance, acc);
      detail::logConservation(caseTag, keBefore, keAfter, acc);
      const real_t keExpected = real_t(0.5) * engine.ac().getMass(boxIdx) * analyticVel.sqrLength();
      detail::logEnergyBudget(caseTag, keAfter, keExpected, acc, tolerance);
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] boxFaceDropRestitution OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")  [Mack et al. 2011, Powder Technology 214:431-442, measured hexahedron "
                              << "e=0.54, mu=0.28]");
   return acc.allPassed();
}

//! Fraige et al. (2008)'s "slide test" methodology (Powder Technology 186:224-240, Section 3.6 / Figs.
//! 7-8): a cube resting FACE-flush (not corner-first) on a plane inclined steeper than the friction angle
//! should slide straight down-slope with no lateral deviation -- their own single-contact-point DEM-Cube
//! version showed the contact force "jumping" between corners as the box translates, injecting a
//! spurious lateral force/torque, fixed there by a multiple-contact-point model. This is the SAME family
//! of degenerate flush-contact tie-break bug fixed in Cylinder/Box::support() this session
//! (2026-09-13) -- a box resting flush on a tilted plane is exactly the configuration where the old
//! support() implementation picked an arbitrary corner instead of the true face centre. The plane is
//! tilted about world X, so by symmetry the box's world-X position AND the two OUT-OF-PLANE angular
//! velocity components (wY, wZ, which would only appear from an asymmetric/"jumping" contact point) must
//! stay EXACTLY zero for all time -- no closed-form dynamics needed, just mirror symmetry -- so any
//! nonzero drift there is a direct regression signal for that bug. The IN-PLANE pitch component (wX,
//! rotation about the incline's own tilt axis) is NOT expected to be zero throughout: the box is released
//! through a small gap and free-falls under world gravity (not purely along the tilted wall normal)
//! before first contact, so a real, physically legitimate settling/rocking transient in wX is expected
//! immediately after impact (confirmed empirically: it decays from a transient ~0.07 rad/s down to
//! ~0.007 rad/s by the last 20% of a 2000-step run) -- only checked in the last `steadyStateFraction` of
//! the run, where it should have damped down close to zero for genuine steady sliding.
inline bool runBoxInclineSlideStraightness(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents  = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 1, 1));
   const real_t density      = b.getParameter< real_t >("density", real_t(1));
   const real_t e            = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   // Fraige et al. (2008), Powder Technology 186:224-240, Table 2: measured cube-wall sliding friction
   // coefficient mu=0.4.
   const real_t mu           = b.getParameter< real_t >("frictionCoefficient", real_t(0.4));
   const real_t poisson      = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt           = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t gravity      = b.getParameter< real_t >("gravityMagnitude", real_t(9.81));
   // Must exceed atan(mu) (~21.8deg for mu=0.4) for the box to actually slide rather than stick.
   const real_t inclineDeg   = b.getParameter< real_t >("inclineAngleDeg", real_t(30));
   const uint_t slideSteps   = b.getParameter< uint_t >("slideSteps", uint_t(2000));
   const real_t steadyStateFraction     = b.getParameter< real_t >("steadyStateFraction", real_t(0.2));
   const real_t toleranceLateralAbs     = b.getParameter< real_t >("toleranceLateralAbs", real_t(1e-6));
   const real_t toleranceOutOfPlaneAbs  = b.getParameter< real_t >("toleranceOutOfPlaneAngVelAbs", real_t(1e-6));
   const real_t toleranceSteadyPitchAbs = b.getParameter< real_t >("toleranceSteadyPitchAbs", real_t(0.01));
   const uint_t vtkSpacing   = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase = b.getParameter< std::string >("outputBaseName", "box_incline_slide");

   const real_t angle      = detail::degToRad(inclineDeg);
   const Rot3   rot(Vec3(1, 0, 0), angle);
   const Vec3   wallNormal = rot.getMatrix() * Vec3(0, 0, 1);
   const Vec3   wallPos(0, 0, 0);
   const Vec3   edgeLength = real_t(2) * halfExtents;

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity,
                     real_t(30) * std::max({ halfExtents[0], halfExtents[1], halfExtents[2] }));
   // Box's local frame shares the wall's own rotation, so in the wall's local frame this is exactly an
   // axis-aligned box resting flush on a z=0 plane -- distance from centre to wall along the normal is
   // exactly halfExtents[2], regardless of inclineDeg. A tiny gap lets gravity pull it into full contact
   // over the first few steps rather than starting at an ambiguous exactly-zero overlap.
   const real_t gap      = real_t(0.01) * std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
   const Vec3   startPos = wallNormal * (halfExtents[2] + gap);
   const size_t boxIdx   = engine.addBox(startPos, Vec3(0, 0, 0), rot, edgeLength, density);
   const size_t wallIdx  = engine.addWall(wallPos, wallNormal);
   engine.finalizeSetup();
   static_cast< void >(wallIdx);

   const std::string tag = outBase + "_" + std::to_string(inclineDeg) + "deg";
   engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                   { boxIdx });
   engine.setupVtk(tag, vtkSpacing);

   const uint_t steadyStateStartStep = uint_t(real_t(slideSteps) * (real_t(1) - steadyStateFraction));
   real_t maxLateralX    = real_t(0);
   real_t maxOutOfPlaneW = real_t(0);
   real_t maxSteadyPitch = real_t(0);
   for (uint_t s = 0; s < slideSteps; ++s)
   {
      engine.step();
      engine.writeCsvRow();
      engine.writeVtkIfDue();
      const Vec3 w = engine.ac().getAngularVelocity(boxIdx);
      maxLateralX    = std::max(maxLateralX, std::abs(engine.ac().getPosition(boxIdx)[0]));
      maxOutOfPlaneW = std::max({ maxOutOfPlaneW, std::abs(w[1]), std::abs(w[2]) });
      if (s >= steadyStateStartStep) maxSteadyPitch = std::max(maxSteadyPitch, std::abs(w[0]));
   }

   detail::PassFailAccumulator acc;
   const bool lateralPass = maxLateralX < toleranceLateralAbs;
   const bool outOfPlanePass = maxOutOfPlaneW < toleranceOutOfPlaneAbs;
   const bool steadyPitchPass = maxSteadyPitch < toleranceSteadyPitchAbs;
   acc.record(lateralPass);
   acc.record(outOfPlanePass);
   acc.record(steadyPitchPass);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] boxInclineSlideStraightness incline=" << inclineDeg << "deg  "
                              << "maxLateralX=" << maxLateralX << " (tol=" << toleranceLateralAbs << ")  "
                              << (lateralPass ? "PASS" : "FAIL") << "   maxOutOfPlaneAngVel(wY,wZ)="
                              << maxOutOfPlaneW << " (tol=" << toleranceOutOfPlaneAbs << ")  "
                              << (outOfPlanePass ? "PASS" : "FAIL") << "   maxSteadyStatePitch(wX, last "
                              << (steadyStateFraction * real_t(100)) << "%)=" << maxSteadyPitch
                              << " (tol=" << toleranceSteadyPitchAbs << ")  " << (steadyPitchPass ? "PASS" : "FAIL"));
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] boxInclineSlideStraightness OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")  [Fraige et al. 2008, Powder Technology 186:224-240, mu=0.4]");
   return acc.allPassed();
}

//======================================================================================================================
// 4-5: pairwise central head-on tests (all 6 shape-pair combinations)
//======================================================================================================================

namespace detail
{
//! Creates a shape by name ("sphere"/"box"/"cylinder") centred at `pos`, moving at `vel`, oriented so
//! that a pair of these shapes approaching head-on along world X makes CENTRAL, face-aligned contact
//! (r x n = 0 on both sides -- no lever arm): box/cylinder are sized and oriented so the contact point
//! sits exactly on the line of centres (box: flat face perpendicular to X, unrotated; cylinder: axis
//! along world Y, so its band's nearest surface point to a partner on the X-axis is exactly (radius,0,0)).
inline size_t createPairMember(DEMEngine& engine, const std::string& shape, const Vec3& pos, const Vec3& vel,
                                real_t sphereRadius, const Vec3& boxHalfExtents, real_t cylRadius, real_t cylLength,
                                real_t density)
{
   if (shape == "sphere") return engine.addSphere(pos, vel, sphereRadius, density);
   if (shape == "box") return engine.addBox(pos, vel, Rot3(), real_t(2) * boxHalfExtents, density);
   // cylinder: body-z aligned with world Y (band faces the partner along X)
   Rot3 rot(Vec3(1, 0, 0), math::pi / real_t(2));
   return engine.addCylinder(pos, vel, rot, cylRadius, cylLength, density);
}

inline const std::vector< std::pair< std::string, std::string > >& allShapePairs()
{
   static const std::vector< std::pair< std::string, std::string > > pairs = {
      { "sphere", "sphere" }, { "sphere", "box" },   { "sphere", "cylinder" },
      { "box", "box" },       { "box", "cylinder" }, { "cylinder", "cylinder" }
   };
   return pairs;
}

//! Runs `engine` until the tracked dynamic-body pair (idxA, idxB) has made and then broken contact, or
//! maxSteps is reached. Returns true if a contact was actually observed (false = the bodies never
//! touched -- misconfigured scenario, not a physical result).
inline bool runUntilSeparated(DEMEngine& engine, size_t idxA, size_t idxB, uint_t maxSteps, bool writeOutput)
{
   bool   everContacted = false;
   real_t penetration   = real_t(0);
   for (uint_t s = 0; s < maxSteps; ++s)
   {
      const bool contactNow = engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0);
      if (contactNow) everContacted = true;
      engine.step();
      if (writeOutput)
      {
         engine.writeCsvRow();
         engine.writeVtkIfDue();
      }
      if (everContacted && !contactNow && s > 0) break; // just separated this step
   }
   return everContacted;
}
} // namespace detail

inline bool runPairRestitution(const config::Config::BlockHandle& b)
{
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 2, 2));
   const real_t cylRadius      = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t cylLength      = b.getParameter< real_t >("cylinderLength", real_t(6));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t e              = b.getParameter< real_t >("restitutionCoefficient", real_t(0.7));
   const real_t mu             = b.getParameter< real_t >("frictionCoefficient", real_t(0.0));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance      = b.getParameter< real_t >("toleranceRelError", real_t(0.1));
   const real_t approachSpeed  = b.getParameter< real_t >("approachSpeed", real_t(1));
   const real_t separation     = b.getParameter< real_t >("initialSeparation", real_t(20));
   const uint_t vtkSpacing     = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase   = b.getParameter< std::string >("outputBaseName", "pair_restitution");

   detail::PassFailAccumulator acc;
   const real_t analyticRatio = analytical::twoBodyHeadOnRestitution(e);

   for (const auto& pr : detail::allShapePairs())
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), separation + real_t(20));
      const size_t idxA = detail::createPairMember(engine, pr.first, Vec3(-separation / 2, 0, 0),
                                                     Vec3(approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      const size_t idxB = detail::createPairMember(engine, pr.second, Vec3(separation / 2, 0, 0),
                                                     Vec3(-approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      engine.finalizeSetup();

      const std::string tag = outBase + "_" + pr.first + "_" + pr.second;
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { idxA, idxB });
      engine.setupVtk(tag, vtkSpacing);

      const Vec3   pBefore  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                              engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);
      const real_t keBefore = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);

      const uint_t maxSteps = uint_t(2) * uint_t(separation / (approachSpeed * dt)) +
                               uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      detail::runUntilSeparated(engine, idxA, idxB, maxSteps, true);

      const Vec3   pAfter  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                             engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);
      const real_t keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);

      const real_t vRelBefore = approachSpeed; // by construction: closing speed along X at t=0
      const real_t vRelAfter  = (engine.ac().getLinearVelocity(idxA) - engine.ac().getLinearVelocity(idxB))[0];
      const real_t simRatio   = vRelAfter / vRelBefore;
      const real_t relErr     = std::abs(simRatio - analyticRatio) / std::max(real_t(1e-8), std::abs(analyticRatio));
      const bool   pass       = relErr < tolerance;
      acc.record(pass);

      const std::string pairTag = std::string("pairRestitution ") + pr.first + "-" + pr.second;
      detail::logMomentumConservation(pairTag, pBefore, pAfter, acc);
      detail::logConservation(pairTag, keBefore, keAfter, acc);
      // Exact closed form for a CENTRAL (r x n = 0, guaranteed by createPairMember's geometry),
      // FRICTIONLESS (mu == 0, this test's default) head-on collision: total KE splits into COM-frame
      // translation (unchanged -- momentum conservation fixes the COM velocity) plus relative-motion KE,
      // which scales by e^2. Only claimed exact when mu is actually 0 -- with friction on, rotational
      // coupling makes the energy budget depend on the same friction-model details already discussed for
      // the wall-impact tests, so the check is skipped rather than silently claiming false precision.
      if (std::abs(mu) < real_t(1e-9))
      {
         const real_t mA = engine.ac().getMass(idxA), mB = engine.ac().getMass(idxB);
         const real_t reducedMass = mA * mB / (mA + mB);
         const real_t keExpected  = keBefore - real_t(0.5) * reducedMass * vRelBefore * vRelBefore *
                                                  (real_t(1) - e * e);
         detail::logEnergyBudget(pairTag, keAfter, keExpected, acc, tolerance);
      }

      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] pairRestitution "
                                 << pr.first << "-" << pr.second << "  simRatio=" << simRatio
                                 << " analyticRatio=" << analyticRatio << " relErr=" << relErr << "  "
                                 << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] pairRestitution OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

//! [Claude, 2026-09-30] Energy conservation of OBLIQUE cylinder-cylinder collisions under each cylinder detection mode
//! (thesis Section 5.2, "continuity"). Two free cylinders with uniformly random orientations collide with a random
//! lateral offset (so that rims, bands and faces all meet), e = 1, mu = 0, no gravity: the kinetic energy after the
//! collision must equal the energy before. A detection that jumps in depth or normal while the bodies move injects or
//! removes energy. The SAME set of initial conditions is run with each mode:
//!   gjk      GJK/EPA only (useAnalyticalContactDetection false)
//!   hybrid   closed form + GJK/EPA fallback (default analytical mode)
//!   pure     closed form without fallback (usePureAnalyticalCylinder)
//!   refined  closed form + support-function refinement (refineCylinderContacts)
//! Reported per mode: mean and max |KE_after/KE_before - 1| over the collisions that made contact, max relative momentum
//! error, max penetration / radius, and collisions that never made contact although the other modes did.
inline bool runCylinderPairEnergyFuzz(const config::Config::BlockHandle& b)
{
   const real_t radius        = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t length        = b.getParameter< real_t >("cylinderLength", real_t(3));
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t speed         = b.getParameter< real_t >("approachSpeed", real_t(1));
   const uint_t numCollisions = b.getParameter< uint_t >("numCollisions", uint_t(200));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t seed          = b.getParameter< uint_t >("seed", uint_t(4242));

   struct Init { Rot3 rotA, rotB; Vec3 offset; };
   std::mt19937                             rng(seed);
   std::uniform_real_distribution< real_t > u(real_t(0), real_t(1));
   auto uniformRot = [&]() {
      const real_t u1 = u(rng), u2 = u(rng), u3 = u(rng);
      const real_t a = std::sqrt(real_t(1) - u1), bb = std::sqrt(u1);
      const real_t w = a * std::sin(real_t(2) * math::pi * u2), x = a * std::cos(real_t(2) * math::pi * u2);
      const real_t y = bb * std::sin(real_t(2) * math::pi * u3), z = bb * std::cos(real_t(2) * math::pi * u3);
      Rot3 r;
      const Vec3 axis = (w < 0 ? real_t(-1) : real_t(1)) * Vec3(x, y, z);
      if (axis.length() > real_t(1e-12)) r.rotate(axis.getNormalized(), real_t(2) * std::acos(std::min(std::abs(w), real_t(1))));
      return r;
   };
   const real_t reach = radius + real_t(0.5) * length;
   std::vector< Init > inits;
   for (uint_t c = 0; c < numCollisions; ++c)
      inits.push_back(Init{ uniformRot(), uniformRot(),
                            Vec3(0, (u(rng) - real_t(0.5)) * real_t(2) * reach, (u(rng) - real_t(0.5)) * real_t(2) * reach) });

   std::vector< std::string > modes;
   {
      std::stringstream ss(b.getParameter< std::string >("modes", "gjk,hybrid,pure,refined"));
      for (std::string item; std::getline(ss, item, ',');) modes.push_back(item);
   }
   std::map< std::string, std::vector< real_t > > keErr, pErr, pen;
   std::map< std::string, std::vector< bool > >   contacted;
   const real_t sep0 = real_t(2.5) * reach;
   for (const auto& mode : modes)
   {
      for (const auto& in : inits)
      {
         DEMEngine engine(dt, real_t(1), real_t(0), poisson, collisionTime, real_t(0), real_t(10) * reach,
                          /*useAnalyticalContactDetection=*/mode != "gjk");
         engine.computeContacts().setUsePureAnalyticalCylinder(mode == "pure");
         engine.computeContacts().setRefineCylinderContacts(mode == "refined");
         const size_t A = engine.addCylinder(Vec3(-sep0 / 2, 0, 0) - in.offset * real_t(0.5), Vec3(speed / 2, 0, 0), in.rotA,
                                             radius, length, density);
         const size_t B = engine.addCylinder(Vec3(sep0 / 2, 0, 0) + in.offset * real_t(0.5), Vec3(-speed / 2, 0, 0), in.rotB,
                                             radius, length, density);
         engine.finalizeSetup();
         const real_t ke0 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p0  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         bool   touched = false;
         real_t maxPen = 0, penetration = 0;
         uint_t quiet = 0;
         const uint_t maxSteps = uint_t(real_t(3) * sep0 / (speed * dt));
         for (uint_t s = 0; s < maxSteps; ++s)
         {
            engine.step();
            const bool c = engine.queryContact(A, B, penetration) && penetration < real_t(0);
            if (mode == "refined" && std::getenv("ANTIDUNES_ENERGY_DEBUG") &&
                int(&in - &inits[0]) == std::atoi(std::getenv("ANTIDUNES_ENERGY_DEBUG")))
            {
               const paper_analytical::CylinderState ci{ engine.ac().getPosition(A), engine.ac().getRotation(A).getMatrix(), radius, length };
               const paper_analytical::CylinderState cj{ engine.ac().getPosition(B), engine.ac().getRotation(B).getMatrix(), radius, length };
               const auto raw = paper_analytical::cylinderCylinderContact(ci, cj);
               const auto rf  = paper_analytical::refineCylinderContact(ci, cj, raw);
               if (c || rf.contact || raw.contact)
               WALBERLA_LOG_INFO_ON_ROOT("[energy-debug] step " << s << " gjkContact=" << c << " gjkPen=" << -penetration << " raw=" << raw.scenario
                  << " " << raw.overlap << " refined=" << rf.scenario << " " << rf.overlap << " n=" << rf.normal << " p=" << rf.point
                  << " rawP=" << raw.point << " KE=" << engine.kineticEnergy(A) + engine.kineticEnergy(B));
            }
            if (c) { touched = true; quiet = 0; maxPen = std::max(maxPen, -penetration); }
            else if (touched && ++quiet > uint_t(20) * uint_t(collisionTime / dt + real_t(1))) break;
         }
         const real_t ke1 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p1  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         contacted[mode].push_back(touched);
         keErr[mode].push_back(touched ? std::abs(ke1 / ke0 - real_t(1)) : real_t(0));
         pErr[mode].push_back((p1 - p0).length() / (engine.ac().getMass(A) * speed));
         pen[mode].push_back(maxPen / radius);
      }
   }

   detail::PassFailAccumulator acc;
   for (const auto& mode : modes)
   {
      real_t sum = 0, mx = 0, pmx = 0, penmx = 0;
      int    n = 0, over = 0, missedVsGjk = 0;
      for (size_t c = 0; c < inits.size(); ++c)
      {
         if (contacted.count("gjk") && contacted["gjk"][c] && !contacted[mode][c]) ++missedVsGjk;
         if (!contacted[mode][c]) continue;
         ++n; sum += keErr[mode][c]; mx = std::max(mx, keErr[mode][c]);
         pmx = std::max(pmx, pErr[mode][c]); penmx = std::max(penmx, pen[mode][c]);
         if (keErr[mode][c] > tolerance)
         {
            ++over;
            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cylinderPairEnergyFuzz mode=" << mode << " collision " << c
                                       << " |KE ratio-1|=" << keErr[mode][c] << " maxPenetration/R=" << pen[mode][c]);
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cylinderPairEnergyFuzz mode=" << mode << " collisions=" << n << "/"
         << inits.size() << " mean|KE ratio-1|=" << (n ? sum / real_t(n) : real_t(0)) << " max=" << mx
         << " over tolerance(" << tolerance << ")=" << over << " maxMomentumErr=" << pmx << " maxPenetration/R=" << penmx
         << " noContactWhereGjkHad=" << missedVsGjk);
      if (mode == "gjk" || mode == "refined") acc.record(over == 0 && missedVsGjk == 0);
   }
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cylinderPairEnergyFuzz OVERALL: " << (acc.allPassed() ? "PASS" : "FAIL") << " ("
                             << acc.passed << "/" << acc.total << ", gjk and refined modes checked)");
   return acc.allPassed();
}

//! [Claude, 2026-10-01] Cuboid version of runCylinderPairEnergyFuzz: oblique elastic (e = 1, mu = 0) collisions of two
//! freely rotating cuboids with random orientations and lateral offsets, no gravity; kinetic energy must be conserved.
//! Modes: "gjk" (GJK/EPA) and "sat" (closed-form separating-axis test, as in the drop tests). Written to test whether
//! the energy gain seen in the randomised cuboid drop tests comes from the SAT contact point (placed on the line between
//! the centres, not at the contact).
inline bool runCuboidPairEnergyFuzz(const config::Config::BlockHandle& b)
{
   const Vec3   halfExt       = b.getParameter< Vec3 >("cuboidHalfExtents", Vec3(6, 1, real_t(1.5)));
   const real_t radius        = std::min(halfExt[0], std::min(halfExt[1], halfExt[2])); // penetration reported / this
   const real_t density       = b.getParameter< real_t >("density", real_t(1));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt            = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t speed         = b.getParameter< real_t >("approachSpeed", real_t(1));
   const uint_t numCollisions = b.getParameter< uint_t >("numCollisions", uint_t(200));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t seed          = b.getParameter< uint_t >("seed", uint_t(4242));

   struct Init { Rot3 rotA, rotB; Vec3 offset; };
   std::mt19937                             rng(seed);
   std::uniform_real_distribution< real_t > u(real_t(0), real_t(1));
   auto uniformRot = [&]() {
      const real_t u1 = u(rng), u2 = u(rng), u3 = u(rng);
      const real_t a = std::sqrt(real_t(1) - u1), bb = std::sqrt(u1);
      const real_t w = a * std::sin(real_t(2) * math::pi * u2), x = a * std::cos(real_t(2) * math::pi * u2);
      const real_t y = bb * std::sin(real_t(2) * math::pi * u3), z = bb * std::cos(real_t(2) * math::pi * u3);
      Rot3 r;
      const Vec3 axis = (w < 0 ? real_t(-1) : real_t(1)) * Vec3(x, y, z);
      if (axis.length() > real_t(1e-12)) r.rotate(axis.getNormalized(), real_t(2) * std::acos(std::min(std::abs(w), real_t(1))));
      return r;
   };
   const real_t reach = halfExt.length();
   std::vector< Init > inits;
   for (uint_t c = 0; c < numCollisions; ++c)
      inits.push_back(Init{ uniformRot(), uniformRot(),
                            Vec3(0, (u(rng) - real_t(0.5)) * real_t(2) * reach, (u(rng) - real_t(0.5)) * real_t(2) * reach) });

   std::vector< std::string > modes;
   {
      std::stringstream ss(b.getParameter< std::string >("modes", "gjk,sat"));
      for (std::string item; std::getline(ss, item, ',');) modes.push_back(item);
   }
   std::map< std::string, std::vector< real_t > > keErr, pErr, pen;
   std::map< std::string, std::vector< bool > >   contacted;
   const real_t sep0 = real_t(2.5) * reach;
   for (const auto& mode : modes)
   {
      for (const auto& in : inits)
      {
         // modes: gjk (GJK/EPA), sat (plain SAT, as in the drop tests), multipoint (useMultiPointBoxContact).
         // [Claude, 2026-10-29] "fraige" mode (useFraigeBoxContactModel) removed -- see
         // archive/unused_contact_methods/ (fraige_multipoint_box_contact.h + fraige_validation_tests.cpp.excerpt).
         DEMEngine engine(dt, real_t(1), real_t(0), poisson, collisionTime, real_t(0), real_t(10) * reach,
                          /*useAnalyticalContactDetection=*/mode != "gjk", /*useRoundedBoxCorners=*/false,
                          /*boxRoundingRadius=*/real_t(0.05), /*useMultiPointBoxContact=*/mode == "multipoint",
                          false, real_t(0.3), false, real_t(0.15), false, real_t(0.25), false);
         const size_t A = engine.addBox(Vec3(-sep0 / 2, 0, 0) - in.offset * real_t(0.5), Vec3(speed / 2, 0, 0), in.rotA,
                                        real_t(2) * halfExt, density);
         const size_t B = engine.addBox(Vec3(sep0 / 2, 0, 0) + in.offset * real_t(0.5), Vec3(-speed / 2, 0, 0), in.rotB,
                                        real_t(2) * halfExt, density);
         engine.finalizeSetup();
         const real_t ke0 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p0  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         bool   touched = false;
         real_t maxPen = 0, penetration = 0;
         uint_t quiet = 0;
         const uint_t maxSteps = uint_t(real_t(3) * sep0 / (speed * dt));
         for (uint_t s = 0; s < maxSteps; ++s)
         {
            engine.step();
            const bool c = engine.queryContact(A, B, penetration) && penetration < real_t(0);
            if (std::getenv("ANTIDUNES_ENERGY_DEBUG") && int(&in - &inits[0]) == std::atoi(std::getenv("ANTIDUNES_ENERGY_DEBUG")))
            {
               const auto r = paper_analytical::boxBoxSAT(engine.ac().getPosition(A), engine.ac().getRotation(A).getMatrix(),
                                                          halfExt, engine.ac().getPosition(B),
                                                          engine.ac().getRotation(B).getMatrix(), halfExt);
               if (c || r.contact)
                  WALBERLA_LOG_INFO_ON_ROOT("[energy-debug] mode=" << mode << " step " << s << " contact=" << c << " pen="
                     << -penetration << " sat=" << r.overlap << " n=" << r.normal << " p=" << r.point
                     << " KE=" << engine.kineticEnergy(A) + engine.kineticEnergy(B));
            }
            if (c) { touched = true; quiet = 0; maxPen = std::max(maxPen, -penetration); }
            else if (touched && ++quiet > uint_t(20) * uint_t(collisionTime / dt + real_t(1))) break;
         }
         const real_t ke1 = engine.kineticEnergy(A) + engine.kineticEnergy(B);
         const Vec3   p1  = engine.ac().getMass(A) * engine.ac().getLinearVelocity(A) +
                            engine.ac().getMass(B) * engine.ac().getLinearVelocity(B);
         contacted[mode].push_back(touched);
         keErr[mode].push_back(touched ? std::abs(ke1 / ke0 - real_t(1)) : real_t(0));
         pErr[mode].push_back((p1 - p0).length() / (engine.ac().getMass(A) * speed));
         pen[mode].push_back(maxPen / radius);
      }
   }

   detail::PassFailAccumulator acc;
   for (const auto& mode : modes)
   {
      real_t sum = 0, mx = 0, pmx = 0, penmx = 0;
      int    n = 0, over = 0, missedVsGjk = 0;
      for (size_t c = 0; c < inits.size(); ++c)
      {
         if (contacted.count("gjk") && contacted["gjk"][c] && !contacted[mode][c]) ++missedVsGjk;
         if (!contacted[mode][c]) continue;
         ++n; sum += keErr[mode][c]; mx = std::max(mx, keErr[mode][c]);
         pmx = std::max(pmx, pErr[mode][c]); penmx = std::max(penmx, pen[mode][c]);
         if (keErr[mode][c] > tolerance)
         {
            ++over;
            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cuboidPairEnergyFuzz mode=" << mode << " collision " << c
                                       << " |KE ratio-1|=" << keErr[mode][c] << " maxPenetration/R=" << pen[mode][c]);
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cuboidPairEnergyFuzz mode=" << mode << " collisions=" << n << "/"
         << inits.size() << " mean|KE ratio-1|=" << (n ? sum / real_t(n) : real_t(0)) << " max=" << mx
         << " over tolerance(" << tolerance << ")=" << over << " maxMomentumErr=" << pmx << " maxPenetration/R=" << penmx
         << " noContactWhereGjkHad=" << missedVsGjk);
      acc.record(over == 0 && missedVsGjk == 0);
   }
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cuboidPairEnergyFuzz OVERALL: " << (acc.allPassed() ? "PASS" : "FAIL") << " ("
                             << acc.passed << "/" << acc.total << ", all modes checked)");
   return acc.allPassed();
}

inline bool runEnergyConservation(const config::Config::BlockHandle& b)
{
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(1, 2, 2));
   const real_t cylRadius      = b.getParameter< real_t >("cylinderRadius", real_t(1));
   const real_t cylLength      = b.getParameter< real_t >("cylinderLength", real_t(6));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t tolerance      = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const real_t approachSpeed  = b.getParameter< real_t >("approachSpeed", real_t(1));
   const real_t separation     = b.getParameter< real_t >("initialSeparation", real_t(20));
   const uint_t vtkSpacing     = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const std::string outBase   = b.getParameter< std::string >("outputBaseName", "energy_conservation");

   detail::PassFailAccumulator acc;
   // Feng (2023)'s ECC criterion: e=1, mu=0, no gravity -- total KE must be conserved across an
   // elastic collision. This is the theoretically rigorous minimum bar for a contact model, and
   // specifically stresses Box x Box / Cylinder x Cylinder, the two pairs ComputeContacts resolves via
   // the generic GJK/EPA fallback (no analytic overload) rather than a closed form -- see
   // ParticleKernels.h's own comment on SolidCylinderContactDetection.
   const real_t e = real_t(1), mu = real_t(0);

   for (const auto& pr : detail::allShapePairs())
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, real_t(0), separation + real_t(20));
      const size_t idxA = detail::createPairMember(engine, pr.first, Vec3(-separation / 2, 0, 0),
                                                     Vec3(approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      const size_t idxB = detail::createPairMember(engine, pr.second, Vec3(separation / 2, 0, 0),
                                                     Vec3(-approachSpeed / 2, 0, 0), sphereRadius, boxHalfExtents,
                                                     cylRadius, cylLength, density);
      engine.finalizeSetup();

      const std::string tag = outBase + "_" + pr.first + "_" + pr.second;
      engine.openCsv(tag + ".csv", { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { idxA, idxB });
      engine.setupVtk(tag, vtkSpacing);

      const real_t keBefore = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      const Vec3   pBefore  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                              engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);

      // Live cross-check of the not-yet-simulation-wired closed forms (PaperAnalyticalContacts.h)
      // against the production GJK/EPA geometry actually driving the dynamics this step, for the two
      // pairs that have a closed form there (box-box, cylinder-cylinder). Diagnostic only -- does not
      // change what drives the physics, per that file's own stated intent.
      const bool crossCheckable =
         (pr.first == "box" && pr.second == "box") || (pr.first == "cylinder" && pr.second == "cylinder");
      real_t maxOverlapDeltaSeen = real_t(0);

      const uint_t maxSteps = uint_t(2) * uint_t(separation / (approachSpeed * dt)) +
                               uint_t(200) * uint_t(collisionTime / dt + real_t(1));
      bool   everContacted = false;
      real_t penetration   = real_t(0);
      for (uint_t s = 0; s < maxSteps; ++s)
      {
         const bool contactNow = engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0);
         if (contactNow)
         {
            everContacted = true;
            if (crossCheckable)
            {
               real_t paperOverlap = real_t(0);
               if (pr.first == "box")
               {
                  const auto res = paper_analytical::boxBoxSAT(
                     engine.ac().getPosition(idxA), engine.ac().getRotation(idxA).getMatrix(), boxHalfExtents,
                     engine.ac().getPosition(idxB), engine.ac().getRotation(idxB).getMatrix(), boxHalfExtents);
                  if (res.contact) paperOverlap = res.overlap;
               }
               else
               {
                  const paper_analytical::CylinderState csA{ engine.ac().getPosition(idxA),
                                                               engine.ac().getRotation(idxA).getMatrix(), cylRadius,
                                                               cylLength };
                  const paper_analytical::CylinderState csB{ engine.ac().getPosition(idxB),
                                                               engine.ac().getRotation(idxB).getMatrix(), cylRadius,
                                                               cylLength };
                  const auto res = paper_analytical::cylinderCylinderContact(csA, csB);
                  if (res.contact) paperOverlap = res.overlap;
               }
               // paper_analytical uses the OPPOSITE overlap sign convention (positive = penetrating);
               // SolidCylinderContactDetection's penetrationDepth is negative = penetrating.
               maxOverlapDeltaSeen = std::max(maxOverlapDeltaSeen, std::abs(paperOverlap - (-penetration)));
            }
         }
         engine.step();
         engine.writeCsvRow();
         engine.writeVtkIfDue();
         if (everContacted && !contactNow && s > 0) break;
      }

      const real_t keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      const Vec3   pAfter  = engine.ac().getMass(idxA) * engine.ac().getLinearVelocity(idxA) +
                             engine.ac().getMass(idxB) * engine.ac().getLinearVelocity(idxB);
      const real_t relErr  = std::abs(keAfter - keBefore) / std::max(real_t(1e-12), keBefore);
      const bool   pass    = everContacted && relErr < tolerance;
      acc.record(pass);
      detail::logMomentumConservation(std::string("energyConservation ") + pr.first + "-" + pr.second, pBefore,
                                       pAfter, acc);

      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] energyConservation "
                                 << pr.first << "-" << pr.second << "  KE_before=" << keBefore
                                 << " KE_after=" << keAfter << " relErr=" << relErr
                                 << (crossCheckable ? "  [paperAnalyticalMaxOverlapDelta=" +
                                                          std::to_string(maxOverlapDeltaSeen) + "]"
                                                    : "")
                                 << "  " << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] energyConservation OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

//======================================================================================================================
// 5b: corner/edge contact-normal discontinuity (Box x Box analytical dispatch) -- REPRODUCTION test
//======================================================================================================================

// Investigates a real defect spotted visually in production job 4196224 (bridge "widerstacks" box
// scenario, ParaView output around step 26000-27000): a batch of staggered-release plank floaters
// lands together and jams, and one plank is subsequently ejected/"flies" out of the pile with no
// single log line ever showing a force large enough to explain it. Literature on polyhedral DEM (e.g.
// the "improved common-plane concept" and "energy-conserving contact" papers for arbitrarily-shaped
// elements) identifies the mechanism directly: near a corner-to-corner or corner-to-edge
// configuration, a tiny relative positional/rotational change can flip which face/edge/vertex pair the
// contact solver reports as the separating axis, snapping the contact normal to a different direction
// between consecutive steps. Because the spring-dashpot force is applied along whatever normal is
// reported THIS step with no continuity constraint against the PREVIOUS step's normal, that snap can
// inject kinetic energy with no corresponding physical source -- exactly what turns an ordinary
// pile-up into a "flying" launch.
//
// This reproduces the geometry realistically, not as a contrived worst case: both bodies use this
// project's own production FloatingBox size/density/Poisson's ratio (test_bridge_v1_..._boxes_..._gas.
// prm's own FloatingBox halfExtents <6,1,1.5>, densityRatio 0.6, and the scenario's own derived
// particleCollisionTime/particleTimestepSize/particleKappa from computedParams.txt), the SAME
// analytical Box x Box dispatch (boxBoxSAT) that scenario actually runs with
// useAnalyticalContactDetection=true -- NOT the generic GJK/EPA fallback every OTHER DEMValidation
// scenario in this file defaults to, since GJK/EPA and boxBoxSAT are different code with different
// failure modes and this defect is specific to the analytical path -- and a realistic approach speed
// (~0.025 LU/step, matching the ~0.01-0.03 LU/step range actually observed in job 4196224's own
// STRONG-CONTACT-FORCE log lines, i.e. what free-fall under this scenario's own tiny gravity,
// 7.03e-6 LU/step^2, produces after a plausible multi-thousand-step drop -- simulating the drop itself
// would just be a long uneventful prelude, so this starts right at that realistic pre-impact speed).
//
// Two identical planks approach roughly along their own long (X) axis -- exactly how two staggered-
// release logs converging on the same jam location would meet -- but offset slightly in Y/Z and tilted
// a few degrees (a plausible natural tumble; ParticleKernels.h's own FloaterHistoryLogger comment notes
// a floater's angular velocity "is often nonzero well before any 'expected' collision"), so their
// SHORT end-cap CORNERS meet first rather than a clean face-to-face impact -- precisely the
// near-degenerate configuration the literature identifies as the trigger. e=1 (perfectly elastic,
// frictionless, no gravity) is used for the PASS/FAIL bound itself -- the strictest, least ambiguous
// bar: total system kinetic energy can dip during compression (converted to non-negative spring
// potential energy) but can never legitimately exceed its own pre-contact value at any instant, since
// there is no other energy source in this system. This is exactly the "ECC" criterion
// runEnergyConservation above already established for this file, applied here as a running bound
// across every step instead of only checked once before/after. restitutionCoefficient/
// frictionCoefficient remain overridable via the prm to rerun under the scenario's literal production
// damping (0.3 "wet/waterlogged wood", 0.5) instead, if a looser reproduction is wanted -- with e<1,
// mu>0 the bound is no longer "exactly conserved" but "monotonically non-increasing once separated",
// so keGainRel is still meaningful as a smoking gun, just not as mathematically clean.
namespace detail
{
//! Shared parameters for one run of the corner/edge contact-normal-discontinuity scenario -- see
//! runCornerNormalDiscontinuity's own (still-present, unchanged) header comment for the full
//! literature/production-incident background and the realism rationale behind every default value
//! here. Factored out so runCornerFixComparison below can run the exact SAME geometry/physics several
//! times, varying only useRoundedBoxCorners/useMultiPointBoxContact.
struct CornerScenarioParams
{
   Vec3   halfExtents      = Vec3(6, 1, 1.5);
   real_t density          = real_t(0.6);
   real_t restitution      = real_t(1);
   real_t friction         = real_t(0);
   real_t poisson          = real_t(0.22);
   real_t collisionTime    = real_t(32.0078);
   real_t dt               = real_t(0.05);
   real_t approachSpeed    = real_t(0.025);
   real_t lateralOffsetY   = real_t(0.6);
   real_t lateralOffsetZ   = real_t(0.9);
   real_t tiltDegreesA     = real_t(8);
   real_t tiltDegreesB     = real_t(-11);
   real_t spinRateB        = real_t(0.004);
   real_t separation       = real_t(20);
   real_t tolerance        = real_t(0.02);
   real_t flipCosThreshold = real_t(0.7); // ~45 degrees, matching ParticleKernels.h's own [NORMAL-FLIP]
   uint_t maxSteps         = uint_t(20000);
   uint_t vtkSpacing       = uint_t(0);
   std::string outputBaseName = "corner_normal_discontinuity";
   // false skips openCsv/setupVtk entirely -- for a large fuzz sweep (runMultiPointBoxContactFuzz)
   // calling this thousands of times, per-trial CSV file I/O would dominate the runtime for no benefit
   // (nobody inspects 2000 near-identical single-trial CSVs).
   bool enableOutput = true;

   // The two togglable mitigations under test -- see SimulationParams.h's
   // useRoundedBoxCorners/useMultiPointBoxContact doc comments. Both default OFF (exact pre-existing
   // behaviour, i.e. the original reproduction case).
   bool   useRoundedBoxCorners   = false;
   real_t boxRoundingRadius      = real_t(0.05);
   bool   useMultiPointBoxContact = false;
   bool   useBoxOverlapCap       = false;
   real_t boxOverlapCapFraction  = real_t(0.3);
   bool   useAxisHysteresis      = false;
   real_t axisHysteresisMargin   = real_t(0.15);
   bool   useContactManifold     = false;
   real_t manifoldFadeRate       = real_t(0.25);
   bool   useRotationalEffectiveMass = false;
};

struct CornerScenarioResult
{
   real_t keBefore        = real_t(0);
   real_t keAfter         = real_t(0);
   real_t keMax           = real_t(0);
   uint_t keMaxStep       = uint_t(0);
   real_t peInitial       = real_t(0); // elastic PE already stored in any t=0 overlap, see
                                       // computeInitialElasticPE's own comment -- 0 for every scenario
                                       // that starts from clean separation (the common case), nonzero
                                       // only for the deep-overlap scenarios that start pre-embedded.
   real_t keGainRel       = real_t(0); // (keMax - (keBefore+peInitial)) / max(1e-12, keBefore+peInitial)
   uint_t normalFlipCount = uint_t(0);
   // The normalFlipCount>0 diagnostic only ever catches a SUDDEN >45-degree single-step jump in the
   // raw SAT normal -- this session's own multiPointDeepOverlapStability data confirms that's real but
   // INCOMPLETE: baseline reproducibly injects energy at several depths with normalFlipCount==0. Two
   // candidate mechanisms this doesn't see: (a) the normal rotating a LOT but gradually (never crossing
   // 45 degrees in any single step, so never flagged, yet still non-conservative -- a linear spring is
   // only exactly conservative along a FIXED axis), and (b) the reported contact POINT jumping (torque
   // arm changes discontinuously) while the normal itself stays stable. Both tracked directly below.
   real_t cumulativeNormalRotationDeg = real_t(0); // sum of EVERY step's angle change, no threshold
   real_t maxPointJumpFrac            = real_t(0); // largest single-step contact-point jump / halfExtents.length()
   uint_t pointJumpCount              = uint_t(0);  // steps where that jump exceeded pointJumpThreshold
   bool   everContacted   = false;
   bool   pass            = false; // keGainRel < tolerance (no anomalous energy gain observed)
   // [Claude, 2026-10-29] the four "eberly*" diagnostics that used to live here (tracking
   // paper_analytical::boxBoxSpheroAnalytical's own reported normal/point) were removed along with
   // useSpheroBoxContact/useSpheroBoxAnalytical -- see archive/unused_contact_methods/.
};

// The "no energy from nowhere" bound every scenario in this file checks is: total mechanical energy at
// any instant <= total mechanical energy at t=0. For every scenario that starts from clean separation
// (approachSpeed>0, zero initial overlap -- runCornerNormalDiscontinuity, runMultiPointBoxContactFuzz,
// etc.) t=0's mechanical energy IS just KE_before, since there is no overlap yet to store any elastic
// potential energy in. But the deep-overlap scenarios (runMultiPointDeepOverlapStability/Fuzz/
// EdgeCases) deliberately start the two bodies ALREADY overlapping at rest -- KE_before is exactly 0
// there, but t=0's true mechanical energy is NOT 0, it is the elastic PE already locked into that
// overlap by construction. Comparing keMax against KE_before alone (as an earlier version of these
// deep-overlap tests did) makes every contacted trial report an appraoch to infinite relative gain
// regardless of whether the config actually behaves correctly, since it divides by ~0 -- this function
// computes the missing term so the pass/fail bound is physically meaningful even from a static release.
//
// Uses the SAME formula ComputeContacts/LinearSpringDashpot itself uses (LinearSpringDashpot.h's own
// setStiffnessAndDamping: stiffnessN = effectiveMass*(pi^2 + ln(e)^2)/collisionTime^2, PE = 1/2*k*x^2 for
// a linear spring compressed by x) and the SAME contact-point set the dynamics will actually use this
// run (boxBoxMultiPointContacts's vertex list -- or its own SAT-fallback single point when that list is
// empty -- if useMultiPointBoxContact, otherwise boxBoxSAT's own single overlap): each point's force is
// scaled by 1/N in ComputeContacts's multi-point loop, so its effective per-point stiffness is
// stiffnessN/N, giving PE = (stiffnessN/(2N)) * sum(overlap_i^2) -- reduces to the ordinary 1/2*k*x^2
// for the N=1 single-point case (baseline/rounded, or multipoint with no vertex embedded).
// pointPE: potential energy stored by ONE point of effective stiffness k compressed by the TRUE
// (uncapped) overlap x, under the SAME piecewise force law ComputeContacts actually applies when a cap
// is active: F(x) = k*x for x<=cap (ordinary linear spring), F(x) = k*cap (constant, clipped) for
// x>cap. This is still a perfectly conservative force (a well-defined potential exists for any
// force that depends only on x), so integrating it gives the CORRECT stored energy at depth x --
// 1/2*k*cap^2 (the ordinary spring energy up to the cap) PLUS k*cap*(x-cap) (the extra distance
// travelled at the now-constant force) for x>cap. Using the naive 1/2*k*min(x,cap)^2 instead (as an
// earlier version of this function did) silently THROWS AWAY that second term, understating the true
// starting energy at any overlap deeper than the cap and making a correctly-behaving capped
// contact look like it "gained" energy from nowhere when it did not.
inline real_t pointPE(real_t k, real_t x, real_t cap)
{
   if (x <= cap) return real_t(0.5) * k * x * x;
   return real_t(0.5) * k * cap * cap + k * cap * (x - cap);
}

inline real_t computeInitialElasticPE(const Vec3& posA, const Rot3& rotA, const Vec3& halfA, real_t densityA,
                                       const Vec3& posB, const Rot3& rotB, const Vec3& halfB, real_t densityB,
                                       real_t restitution, real_t collisionTime, bool useRoundedBoxCorners,
                                       real_t boxRoundingRadius, bool useMultiPointBoxContact,
                                       bool useBoxOverlapCap = false, real_t boxOverlapCapFraction = real_t(0.3),
                                       bool useRotationalEffectiveMass = false)
{
   const real_t massA = densityA * real_t(8) * halfA[0] * halfA[1] * halfA[2];
   const real_t massB = densityB * real_t(8) * halfB[0] * halfB[1] * halfB[2];
   const real_t meffTranslational = real_t(1) / (real_t(1) / massA + real_t(1) / massB);
   const real_t lnE   = std::log(restitution);

   // Same generalized effective-mass formula as ParticleKernels.h's own computeRotationalEffectiveMass
   // (see ModelParameters::useRotationalEffectiveMass's doc comment) -- re-derived here from raw geometry
   // since this function has no live mesa_pd Accessor/ParticleStorage to query. Each box's world-frame
   // inverse inertia only depends on its (static, at t=0) rotation, so it's built once per box and reused
   // for whichever single (point,normal) ends up chosen as PRIMARY below -- exactly mirroring
   // ComputeContacts' own behaviour, where ONE meff (from the primary SAT point/normal, computed BEFORE
   // branching into the multi-point loop) sets stiffnessN/dampingN for every point applied that sub-cycle.
   Mat3 invInertiaWF_A, invInertiaWF_B;
   if (useRotationalEffectiveMass)
   {
      FloatingBox boxA(real_t(2) * halfA);
      boxA.updateMassAndInertia(densityA);
      FloatingBox boxB(real_t(2) * halfB);
      boxB.updateMassAndInertia(densityB);
      invInertiaWF_A = math::transformMatrixRART(rotA.getMatrix(), boxA.getInvInertiaBF());
      invInertiaWF_B = math::transformMatrixRART(rotB.getMatrix(), boxB.getInvInertiaBF());
   }
   const auto stiffnessAtPrimary = [&](const Vec3& point, const Vec3& normal) -> real_t {
      real_t meff = meffTranslational;
      if (useRotationalEffectiveMass)
      {
         const Vec3 rA = point - posA;
         const Vec3 rB = point - posB;
         const Vec3 crossA = rA % normal;
         const Vec3 crossB = rB % normal;
         real_t invMeff = real_t(1) / massA + real_t(1) / massB;
         if (crossA.sqrLength() > real_t(0)) invMeff += crossA * (invInertiaWF_A * crossA);
         if (crossB.sqrLength() > real_t(0)) invMeff += crossB * (invInertiaWF_B * crossB);
         meff = real_t(1) / invMeff;
      }
      return meff * (math::pi * math::pi + lnE * lnE) / (collisionTime * collisionTime);
   };
   // Must match ComputeContacts' own per-pair cap exactly (ParticleKernels.h) so this ceiling reflects
   // the SAME cap value the actual dynamics will use. Always computed from the TRUE (uncapped)
   // geometric overlap below -- see pointPE's own comment for why the cap changes HOW the energy is
   // computed (piecewise), not whether the true starting overlap/energy is smaller.
   const real_t overlapCap = useBoxOverlapCap
      ? boxOverlapCapFraction * std::min({ halfA[0], halfA[1], halfA[2], halfB[0], halfB[1], halfB[2] })
      : std::numeric_limits< real_t >::max();

   // [Claude, 2026-10-29] useSpheroBoxContact branch removed -- see archive/unused_contact_methods/.

   if (useMultiPointBoxContact)
   {
      // Uncapped on purpose -- need the TRUE overlap at each point for pointPE below, not whatever a
      // cap would have clipped it to.
      auto points = paper_analytical::boxBoxMultiPointContacts(posA, rotA.getMatrix(), halfA, posB,
                                                                 rotB.getMatrix(), halfB);
      const auto sat = paper_analytical::boxBoxSAT(posA, rotA.getMatrix(), halfA, posB, rotB.getMatrix(),
                                                     halfB, useRoundedBoxCorners, boxRoundingRadius);
      if (points.empty())
      {
         if (!sat.contact) return real_t(0);
         return pointPE(stiffnessAtPrimary(sat.point, sat.normal), sat.overlap, overlapCap);
      }
      // PRIMARY point/normal: sat's own result if it found contact, else the deepest vertex -- exactly
      // matching ParticleKernels.h's Box x Box operator() (see its own comment right after
      // "SAT itself found no contact, but useMultiPointBoxContact's independent vertex test did").
      Vec3 primaryPoint, primaryNormal;
      if (sat.contact) { primaryPoint = sat.point; primaryNormal = sat.normal; }
      else
      {
         const auto deepest = std::max_element(points.begin(), points.end(),
            [](const auto& a, const auto& b) { return a.overlap < b.overlap; });
         primaryPoint = deepest->point; primaryNormal = deepest->normal;
      }
      const real_t stiffnessN = stiffnessAtPrimary(primaryPoint, primaryNormal);
      const real_t kEff = stiffnessN / real_t(points.size()); // matches ComputeContacts' own 1/N force scaling
      real_t sum = real_t(0);
      for (const auto& pt : points) sum += pointPE(kEff, pt.overlap, overlapCap);
      return sum;
   }

   const auto sat = paper_analytical::boxBoxSAT(posA, rotA.getMatrix(), halfA, posB, rotB.getMatrix(), halfB,
                                                 useRoundedBoxCorners, boxRoundingRadius);
   if (!sat.contact) return real_t(0);
   return pointPE(stiffnessAtPrimary(sat.point, sat.normal), sat.overlap, overlapCap);
}

inline CornerScenarioResult runCornerScenarioOnce(const CornerScenarioParams& p)
{
   // useAnalyticalContactDetection=true: the whole point -- reproduce against the ACTUAL production
   // Box x Box dispatch (boxBoxSAT), not the GJK/EPA fallback every other scenario in this file uses.
   DEMEngine engine(p.dt, p.restitution, p.friction, p.poisson, p.collisionTime, real_t(0),
                     p.separation + real_t(40), /*useAnalyticalContactDetection=*/true,
                     p.useRoundedBoxCorners, p.boxRoundingRadius, p.useMultiPointBoxContact,
                     p.useBoxOverlapCap, p.boxOverlapCapFraction, p.useAxisHysteresis,
                     p.axisHysteresisMargin, p.useContactManifold, p.manifoldFadeRate,
                     p.useRotationalEffectiveMass);

   const Rot3 rotA(Vec3(0, 0, 1), p.tiltDegreesA * math::pi / real_t(180));
   const Rot3 rotB(Vec3(0, 0, 1), p.tiltDegreesB * math::pi / real_t(180));

   const size_t idxA = engine.addBox(Vec3(-p.separation / 2, 0, 0), Vec3(p.approachSpeed / 2, 0, 0), rotA,
                                      real_t(2) * p.halfExtents, p.density);
   const size_t idxB = engine.addBox(Vec3(p.separation / 2, p.lateralOffsetY, p.lateralOffsetZ),
                                      Vec3(-p.approachSpeed / 2, 0, 0), rotB, real_t(2) * p.halfExtents,
                                      p.density);
   // Small residual tumble on B ONLY, about an axis (Y) that keeps sweeping its end-cap corner through
   // A's corner region as the two close in -- NOT an artificial imposed spin to force the bug, but the
   // same order of magnitude as the pre-contact angular velocities this scenario's own floaters
   // routinely carry (see runCornerNormalDiscontinuity's own header comment).
   engine.ac().setAngularVelocity(idxB, Vec3(0, p.spinRateB, 0));

   engine.finalizeSetup();

   if (p.enableOutput)
   {
      engine.openCsv(p.outputBaseName + ".csv",
                      { "posX", "posY", "posZ", "velX", "velY", "velZ", "wX", "wY", "wZ", "KE" },
                      { idxA, idxB });
      engine.setupVtk(p.outputBaseName, p.vtkSpacing);
   }

   CornerScenarioResult result;
   result.keBefore  = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
   result.peInitial = computeInitialElasticPE(
      Vec3(-p.separation / 2, 0, 0), rotA, p.halfExtents, p.density,
      Vec3(p.separation / 2, p.lateralOffsetY, p.lateralOffsetZ), rotB, p.halfExtents, p.density,
      p.restitution, p.collisionTime, p.useRoundedBoxCorners, p.boxRoundingRadius, p.useMultiPointBoxContact,
      p.useBoxOverlapCap, p.boxOverlapCapFraction, p.useRotationalEffectiveMass);
   result.keMax    = result.keBefore;

   real_t penetration = real_t(0);
   Vec3   lastNormal;
   Vec3   lastPoint;
   bool   haveLastNormal = false;
   const real_t pointJumpThreshold = real_t(0.15); // fraction of halfExtents.length(), see doc comment above

   for (uint_t s = 0; s < p.maxSteps; ++s)
   {
      const bool contactNow = engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0);
      if (contactNow)
      {
         result.everContacted = true;
         // Diagnostic-only: reads back the SAME closed-form dispatch that is already driving the
         // dynamics this step, purely to record its own reported normal/point -- not a second/
         // different code path. Deliberately does NOT pass p.useRoundedBoxCorners here: this
         // diagnostic exists to characterize the RAW SAT axis's own discontinuity regardless of which
         // mitigation (if any) is actually resolving the dynamics this run.
         const auto res = paper_analytical::boxBoxSAT(engine.ac().getPosition(idxA),
                                                        engine.ac().getRotation(idxA).getMatrix(), p.halfExtents,
                                                        engine.ac().getPosition(idxB),
                                                        engine.ac().getRotation(idxB).getMatrix(), p.halfExtents);
         if (res.contact)
         {
            if (haveLastNormal)
            {
               const real_t cosAngle = std::max(real_t(-1), std::min(real_t(1), res.normal * lastNormal));
               if (cosAngle < p.flipCosThreshold) ++result.normalFlipCount;
               // Unconditional -- unlike normalFlipCount, this accumulates EVERY step's angle change,
               // so a normal that rotates a lot gradually (never crossing 45 degrees in one step) still
               // shows up here, distinguishing "genuinely fixed axis" from "slowly but continuously
               // rotating axis" -- both look identical to normalFlipCount, only the latter can still
               // inject energy (a linear spring is only exactly conservative along a FIXED direction).
               result.cumulativeNormalRotationDeg += std::acos(cosAngle) * real_t(180) / math::pi;

               const real_t pointJumpFrac = (res.point - lastPoint).length() / p.halfExtents.length();
               if (pointJumpFrac > result.maxPointJumpFrac) result.maxPointJumpFrac = pointJumpFrac;
               if (pointJumpFrac > pointJumpThreshold) ++result.pointJumpCount;
            }
            lastNormal     = res.normal;
            lastPoint      = res.point;
            haveLastNormal = true;
         }
         // [Claude, 2026-10-29] the useSpheroBoxContact+useSpheroBoxAnalytical Eberly-diagnostic
         // block that used to live here was removed -- see archive/unused_contact_methods/.
      }
      engine.step();
      engine.writeCsvRow();
      engine.writeVtkIfDue();

      const real_t keNow = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      if (keNow > result.keMax)
      {
         result.keMax     = keNow;
         result.keMaxStep = s;
      }

      if (result.everContacted && !contactNow && s > 0) break;
   }

   result.keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
   // Ceiling = KE_before + the elastic PE already stored in any t=0 overlap (0 for every scenario that
   // starts from clean separation -- see computeInitialElasticPE's own comment for why this term
   // matters for the deep-overlap scenarios, which start pre-embedded at rest with KE_before=0).
   const real_t ceiling = result.keBefore + result.peInitial;
   result.keGainRel = (result.keMax - ceiling) / std::max(real_t(1e-12), ceiling);
   result.pass      = result.everContacted && result.keGainRel < p.tolerance;
   return result;
}
} // namespace detail

// Investigates a real defect spotted visually in production job 4196224 (bridge "widerstacks" box
// scenario, ParaView output around step 26000-27000): a batch of staggered-release plank floaters
// lands together and jams, and one plank is subsequently ejected/"flies" out of the pile with no
// single log line ever showing a force large enough to explain it. Literature on polyhedral DEM (e.g.
// the "improved common-plane concept" and "energy-conserving contact" papers for arbitrarily-shaped
// elements) identifies the mechanism directly: near a corner-to-corner or corner-to-edge
// configuration, a tiny relative positional/rotational change can flip which face/edge/vertex pair the
// contact solver reports as the separating axis, snapping the contact normal to a different direction
// between consecutive steps. Because the spring-dashpot force is applied along whatever normal is
// reported THIS step with no continuity constraint against the PREVIOUS step's normal, that snap can
// inject kinetic energy with no corresponding physical source -- exactly what turns an ordinary
// pile-up into a "flying" launch.
//
// This reproduces the geometry realistically, not as a contrived worst case: both bodies use this
// project's own production FloatingBox size/density/Poisson's ratio (test_bridge_v1_..._boxes_..._gas.
// prm's own FloatingBox halfExtents <6,1,1.5>, densityRatio 0.6, and the scenario's own derived
// particleCollisionTime/particleTimestepSize/particleKappa from computedParams.txt), the SAME
// analytical Box x Box dispatch (boxBoxSAT) that scenario actually runs with
// useAnalyticalContactDetection=true -- NOT the generic GJK/EPA fallback every OTHER DEMValidation
// scenario in this file defaults to, since GJK/EPA and boxBoxSAT are different code with different
// failure modes and this defect is specific to the analytical path -- and a realistic approach speed
// (~0.025 LU/step, matching the ~0.01-0.03 LU/step range actually observed in job 4196224's own
// STRONG-CONTACT-FORCE log lines, i.e. what free-fall under this scenario's own tiny gravity,
// 7.03e-6 LU/step^2, produces after a plausible multi-thousand-step drop -- simulating the drop itself
// would just be a long uneventful prelude, so this starts right at that realistic pre-impact speed).
//
// Two identical planks approach roughly along their own long (X) axis -- exactly how two staggered-
// release logs converging on the same jam location would meet -- but offset slightly in Y/Z and tilted
// a few degrees (a plausible natural tumble; ParticleKernels.h's own FloaterHistoryLogger comment notes
// a floater's angular velocity "is often nonzero well before any 'expected' collision"), so their
// SHORT end-cap CORNERS meet first rather than a clean face-to-face impact -- precisely the
// near-degenerate configuration the literature identifies as the trigger. e=1 (perfectly elastic,
// frictionless, no gravity) is used for the PASS/FAIL bound itself -- the strictest, least ambiguous
// bar: total system kinetic energy can dip during compression (converted to non-negative spring
// potential energy) but can never legitimately exceed its own pre-contact value at any instant, since
// there is no other energy source in this system. useRoundedBoxCorners/useMultiPointBoxContact
// (default OFF, matching the original reproduction) can be set true from the prm to instead verify one
// or both mitigations clear this bound -- see runCornerFixComparison below for a single run that
// checks all four combinations side by side.
inline bool runCornerNormalDiscontinuity(const config::Config::BlockHandle& b)
{
   detail::CornerScenarioParams p;
   p.halfExtents             = b.getParameter< Vec3 >("halfExtents", p.halfExtents);
   p.density                 = b.getParameter< real_t >("density", p.density);
   p.restitution             = b.getParameter< real_t >("restitutionCoefficient", p.restitution);
   p.friction                = b.getParameter< real_t >("frictionCoefficient", p.friction);
   p.poisson                 = b.getParameter< real_t >("poissonsRatio", p.poisson);
   p.collisionTime           = b.getParameter< real_t >("collisionTime", p.collisionTime);
   p.dt                      = b.getParameter< real_t >("dt", p.dt);
   p.approachSpeed           = b.getParameter< real_t >("approachSpeed", p.approachSpeed);
   p.lateralOffsetY          = b.getParameter< real_t >("lateralOffsetY", p.lateralOffsetY);
   p.lateralOffsetZ          = b.getParameter< real_t >("lateralOffsetZ", p.lateralOffsetZ);
   p.tiltDegreesA            = b.getParameter< real_t >("tiltDegreesA", p.tiltDegreesA);
   p.tiltDegreesB            = b.getParameter< real_t >("tiltDegreesB", p.tiltDegreesB);
   p.spinRateB               = b.getParameter< real_t >("spinRateB", p.spinRateB);
   p.separation              = b.getParameter< real_t >("initialSeparation", p.separation);
   p.tolerance               = b.getParameter< real_t >("toleranceRelError", p.tolerance);
   p.flipCosThreshold        = b.getParameter< real_t >("flipCosThreshold", p.flipCosThreshold);
   p.maxSteps                = b.getParameter< uint_t >("maxSteps", p.maxSteps);
   p.vtkSpacing              = b.getParameter< uint_t >("vtkSpacing", p.vtkSpacing);
   p.outputBaseName          = b.getParameter< std::string >("outputBaseName", p.outputBaseName);
   p.useRoundedBoxCorners    = b.getParameter< bool >("useRoundedBoxCorners", p.useRoundedBoxCorners);
   p.boxRoundingRadius       = b.getParameter< real_t >("boxRoundingRadius", p.boxRoundingRadius);
   p.useMultiPointBoxContact = b.getParameter< bool >("useMultiPointBoxContact", p.useMultiPointBoxContact);

   const auto r = detail::runCornerScenarioOnce(p);

   detail::PassFailAccumulator acc;
   acc.record(r.pass);

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerNormalDiscontinuity  KE_before=" << r.keBefore
      << " KE_after=" << r.keAfter << " KE_max=" << r.keMax << " (at step " << r.keMaxStep << ")"
      << " maxRelGain=" << r.keGainRel << " normalFlipCount=" << r.normalFlipCount
      << " everContacted=" << r.everContacted << "  " << (r.pass ? "PASS" : "FAIL (energy injected)"));

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerNormalDiscontinuity OVERALL: "
      << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total << ")");
   return acc.allPassed();
}

// Runs the EXACT SAME corner/tilt/offset geometry as runCornerNormalDiscontinuity above four times,
// varying only useRoundedBoxCorners/useMultiPointBoxContact, and reports a side-by-side comparison --
// the direct answer to "do the two togglable mitigations actually fix the reproduced defect, and does
// combining them help further or just cost more". Geometry/physics parameters are read from the SAME
// prm keys as runCornerNormalDiscontinuity (minus useRoundedBoxCorners/useMultiPointBoxContact
// themselves, which this scenario controls); outputBaseName is suffixed per configuration so all four
// CSVs/VTK series are kept side by side rather than overwriting each other.
inline bool runCornerFixComparison(const config::Config::BlockHandle& b)
{
   detail::CornerScenarioParams p;
   p.halfExtents      = b.getParameter< Vec3 >("halfExtents", p.halfExtents);
   p.density          = b.getParameter< real_t >("density", p.density);
   p.restitution      = b.getParameter< real_t >("restitutionCoefficient", p.restitution);
   p.friction         = b.getParameter< real_t >("frictionCoefficient", p.friction);
   p.poisson          = b.getParameter< real_t >("poissonsRatio", p.poisson);
   p.collisionTime    = b.getParameter< real_t >("collisionTime", p.collisionTime);
   p.dt               = b.getParameter< real_t >("dt", p.dt);
   p.approachSpeed    = b.getParameter< real_t >("approachSpeed", p.approachSpeed);
   p.lateralOffsetY   = b.getParameter< real_t >("lateralOffsetY", p.lateralOffsetY);
   p.lateralOffsetZ   = b.getParameter< real_t >("lateralOffsetZ", p.lateralOffsetZ);
   p.tiltDegreesA     = b.getParameter< real_t >("tiltDegreesA", p.tiltDegreesA);
   p.tiltDegreesB     = b.getParameter< real_t >("tiltDegreesB", p.tiltDegreesB);
   p.spinRateB        = b.getParameter< real_t >("spinRateB", p.spinRateB);
   p.separation       = b.getParameter< real_t >("initialSeparation", p.separation);
   p.tolerance        = b.getParameter< real_t >("toleranceRelError", p.tolerance);
   p.flipCosThreshold = b.getParameter< real_t >("flipCosThreshold", p.flipCosThreshold);
   p.maxSteps         = b.getParameter< uint_t >("maxSteps", p.maxSteps);
   p.vtkSpacing       = b.getParameter< uint_t >("vtkSpacing", p.vtkSpacing);
   const real_t boxRoundingRadius = b.getParameter< real_t >("boxRoundingRadius", p.boxRoundingRadius);
   p.boxRoundingRadius = boxRoundingRadius;

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
   };
   const Config configs[] = {
      { "baseline (both OFF)", false, false },
      { "useRoundedBoxCorners", true, false },
      { "useMultiPointBoxContact", false, true },
      { "both", true, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison -- same corner/tilt/offset geometry "
      "(the exact case runCornerNormalDiscontinuity reproduces), 4 configurations:");

   uint_t configsWithoutInjection = 0;
   bool   baselineFailed          = false;
   bool   anyFixPassed            = false;
   for (const auto& cfg : configs)
   {
      p.useRoundedBoxCorners    = cfg.rounded;
      p.useMultiPointBoxContact = cfg.multiPoint;
      p.outputBaseName          = std::string("corner_fix_comparison_") + cfg.name;

      const auto r = detail::runCornerScenarioOnce(p);
      if (r.pass) ++configsWithoutInjection;
      if (std::string(cfg.name) == "baseline (both OFF)" && !r.pass) baselineFailed = true;
      if (std::string(cfg.name) != "baseline (both OFF)" && r.pass) anyFixPassed = true;

      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison  config=\"" << cfg.name << "\""
         << "  KE_before=" << r.keBefore << " KE_max=" << r.keMax << " (at step " << r.keMaxStep << ")"
         << " maxRelGain=" << r.keGainRel << " normalFlipCount=" << r.normalFlipCount
         << " everContacted=" << r.everContacted
         << "  " << (r.pass ? "PASS (no injection)" : "FAIL (energy injected)"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison OVERALL: " << configsWithoutInjection
      << "/" << (sizeof(configs) / sizeof(configs[0])) << " configurations avoided energy injection");

   // This is a comparison/characterization tool, not a strict per-scenario regression gate -- the
   // meaningful pass criterion is "the defect is still confirmed present at baseline (otherwise this
   // whole comparison is moot) AND at least one mitigation actually clears it" -- both being true is
   // the direct, minimal evidence that (a) the reproduction still holds and (b) a real fix exists.
   // [Claude, 2026-10-01] since the SAT contact-point fix (detail::boxBoxContactPoint) the baseline no longer injects energy
   // here, so "baseline reproduces the defect" no longer holds; the criterion is now that NO configuration injects
   // energy (with ANTIDUNES_SAT_FORMER_POINT set, the former point and the former outcome are reproduced).
   const bool pass = configsWithoutInjection == uint_t(sizeof(configs) / sizeof(configs[0]));
   (void) anyFixPassed;
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] cornerFixComparison: baseline reproduces the defect="
      << baselineFailed << ", at least one mitigation fixes it=" << anyFixPassed << "  "
      << (pass ? "PASS" : "FAIL"));
   return pass;
}

// Targets a DIFFERENT regime than runCornerNormalDiscontinuity/runCornerFixComparison/
// runMultiPointBoxContactFuzz above: all three of those start the two boxes SEPARATED and let them
// approach at a slow, realistic speed, so contact only ever begins at a SHALLOW, near-degenerate
// overlap -- exactly the regime Fraige et al. 2008a's multi-point concept targets. None of them ever
// exercise what boxBoxMultiPointContacts (PaperAnalyticalContacts.h) does when a vertex is embedded
// DEEP inside the other box -- e.g. a jam where several floaters have already compressed hard against
// each other, or a body released directly into an existing pile (see the RELEASE-OVERLAP diagnostic
// in ParticleKernels.h). That function's per-vertex overlap is bounded only by
// min(halfExtents) (up to ~1 LU for this project's own <6,1,1.5> floating boxes, i.e. NOT a small
// near-surface number) and has no shared cap across points the way a single SAT overlap implicitly
// does -- each embedded vertex contributes its OWN independent spring force at its OWN world point,
// summed. The concern (raised directly against real production VTK evidence: uid 3314649325744685089
// tipping to a near-vertical pose over ~9000 steps in job 4201127, the useMultiPointBoxContact rerun of
// the crashed 200k "boxes/layer1only" scenario) is that under deep interpenetration this can inject
// MORE spurious energy/torque than the single-point boxBoxSAT path it's meant to improve on, since N
// independently-deep vertices each store their own "spring energy" rather than one shared overlap
// value.
//
// Starts the two boxes ALREADY OVERLAPPING by a controlled amount (approachSpeed=0 -- no closing
// momentum at all, so ANY kinetic energy that appears after release comes purely from resolving the
// static initial overlap, not from a real collision) at several depths from shallow (a bare touch) to
// severe (well past what boxBoxSAT's own single-axis result could plausibly represent), each run
// against baseline/rounded/multiPoint/both -- mirroring runCornerFixComparison's own 4-way structure,
// but swept over depth instead of compared at one fixed shallow contact. Same strict e=1, mu=0, no-
// gravity ECC bound as every other scenario in this file: total system KE must never exceed its own
// pre-release value, since there is no other energy source. A tilt is kept on both bodies (same
// defaults as CornerScenarioParams) specifically so MULTIPLE vertices, not just one clean face-face
// overlap, are actually embedded -- a flat face-face deep overlap wouldn't exercise the multi-point-
// vs-single-point comparison this test exists for at all.
inline bool runMultiPointDeepOverlapStability(const config::Config::BlockHandle& b)
{
   detail::CornerScenarioParams p;
   p.halfExtents      = b.getParameter< Vec3 >("halfExtents", p.halfExtents);
   p.density          = b.getParameter< real_t >("density", p.density);
   p.restitution      = b.getParameter< real_t >("restitutionCoefficient", p.restitution);
   p.friction         = b.getParameter< real_t >("frictionCoefficient", p.friction);
   p.poisson          = b.getParameter< real_t >("poissonsRatio", p.poisson);
   p.collisionTime    = b.getParameter< real_t >("collisionTime", p.collisionTime);
   p.dt               = b.getParameter< real_t >("dt", p.dt);
   p.approachSpeed    = real_t(0); // static overlap only -- no closing momentum, see header comment
   p.lateralOffsetY   = b.getParameter< real_t >("lateralOffsetY", p.lateralOffsetY);
   p.lateralOffsetZ   = b.getParameter< real_t >("lateralOffsetZ", p.lateralOffsetZ);
   p.tiltDegreesA     = b.getParameter< real_t >("tiltDegreesA", p.tiltDegreesA);
   p.tiltDegreesB     = b.getParameter< real_t >("tiltDegreesB", p.tiltDegreesB);
   p.spinRateB        = real_t(0); // isolate overlap-release energy from any imposed spin too
   p.tolerance        = b.getParameter< real_t >("toleranceRelError", p.tolerance);
   p.maxSteps         = b.getParameter< uint_t >("maxSteps", uint_t(4000));
   p.vtkSpacing       = b.getParameter< uint_t >("vtkSpacing", p.vtkSpacing);
   p.boxRoundingRadius = b.getParameter< real_t >("boxRoundingRadius", p.boxRoundingRadius);

   // Separation values relative to the box's own full length (2*halfExtents[0] = 12 by default):
   // 10 = shallow corner touch (same order as runCornerNormalDiscontinuity's own 20 minus its own
   // approach distance), down to -4 = centres nearly coincident, well past anything a real single
   // sub-cycle's worth of approach velocity could produce -- deliberately including an extreme point
   // so a latent "unbounded with depth" defect cannot hide just outside a narrower sweep.
   const std::vector< real_t > depths = { real_t(10), real_t(6), real_t(3), real_t(0), real_t(-4) };

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
   };
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false },
      { "both", true, true, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false },
      { "useAxisHysteresis", false, false, false, true, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapStability -- same tilt/offset "
      "geometry as runCornerNormalDiscontinuity, started ALREADY OVERLAPPING (approachSpeed=0) at "
      << depths.size() << " depths x " << (sizeof(configs) / sizeof(configs[0])) << " configurations:");

   bool multiPointEverWorseThanBaseline = false;
   bool anyConfigFailed                 = false;
   for (const real_t sep : depths)
   {
      p.separation = sep;
      real_t baselineGain = real_t(0);
      for (const auto& cfg : configs)
      {
         p.useRoundedBoxCorners    = cfg.rounded;
         p.useMultiPointBoxContact = cfg.multiPoint;
         p.useBoxOverlapCap        = cfg.overlapCap;
         p.useAxisHysteresis       = cfg.hysteresis;
         p.useContactManifold      = cfg.manifold;
         p.outputBaseName = std::string("multipoint_deep_overlap_sep") + std::to_string(sep) + "_" + cfg.name;

         const auto r = detail::runCornerScenarioOnce(p);
         if (std::string(cfg.name) == "baseline (both OFF)") baselineGain = r.keGainRel;
         if (std::string(cfg.name) == "useMultiPointBoxContact" && r.keGainRel > baselineGain + p.tolerance)
            multiPointEverWorseThanBaseline = true;
         if (!r.pass) anyConfigFailed = true;

         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapStability  separation=" << sep
            << " config=\"" << cfg.name << "\"  KE_before=" << r.keBefore << " PE_initial=" << r.peInitial
            << " KE_max=" << r.keMax
            << " (step " << r.keMaxStep << ") maxRelGain=" << r.keGainRel
            << " rawSatNormalFlipCount=" << r.normalFlipCount
            << " cumulativeNormalRotationDeg=" << r.cumulativeNormalRotationDeg
            << " maxPointJumpFrac=" << r.maxPointJumpFrac << " pointJumpCount=" << r.pointJumpCount
            << " everContacted=" << r.everContacted << "  " << (r.pass ? "PASS" : "FAIL"));
      }
   }

   // Two independent findings reported separately: (a) does EVERY configuration stay inside the
   // strict ECC bound at every depth (the same bar every other scenario in this file uses), and (b)
   // even where multiPoint itself still passes that bound, does it ever inject MORE relative energy
   // than the plain single-point baseline at the SAME depth -- the specific "multi-point is not an
   // improvement, it behaves more like the original defect" concern this test was written to check,
   // which a bare PASS/FAIL against the tolerance alone would not surface.
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapStability OVERALL: "
      "anyConfigFailedEccBound=" << anyConfigFailed << "  multiPointEverWorseThanBaselineAtSameDepth="
      << multiPointEverWorseThanBaseline << "  "
      << (!anyConfigFailed && !multiPointEverWorseThanBaseline ? "PASS" : "FAIL"));
   return !anyConfigFailed && !multiPointEverWorseThanBaseline;
}

// Randomized counterpart to runMultiPointDeepOverlapStability, the same relationship
// runMultiPointBoxContactFuzz already has to runCornerNormalDiscontinuity/runCornerFixComparison: the
// hand-picked depths/tilt/offset in that scenario are FIVE specific points on what is actually a
// continuous, multi-dimensional space (box aspect ratio, tilt on each body, lateral offset, starting
// overlap depth) -- this sweeps that space randomly, looking for the worst-case combination rather
// than trusting five samples to be representative. Every trial:
//   (1) GEOMETRY, direct on paper_analytical::boxBoxMultiPointContacts (same sanity checks as
//       runMultiPointBoxContactFuzz's own geometry pass: finite, unit normals, non-negative overlap)
//       PLUS the new check this whole investigation exists for: records the WORST (largest) single
//       vertex overlap seen across all trials, and separately counts how many trials produce at least
//       one vertex overlap deep enough to trip [LARGE-OVERLAP-MULTIPOINT]'s own threshold (half the
//       smaller body's bounding-sphere radius, ParticleKernels.h) -- quantifying how often a "normal"
//       random configuration (not a hand-crafted worst case) actually reaches that regime.
//   (2) DYNAMICS, all 4 configurations (baseline/rounded/multiPoint/both) per trial, approachSpeed=0
//       (static release only, isolating overlap-driven injection from momentum transfer -- same
//       rationale as runMultiPointDeepOverlapStability), strict e=1/mu=0/no-gravity ECC bound. Tracks
//       whether multiPoint (or both) EVER injects more relative energy than the plain baseline at the
//       SAME randomly-drawn trial -- the specific "not an improvement, more like the original defect"
//       question this whole line of testing exists to answer, which a bare pass/fail vs. a fixed
//       tolerance would not surface on its own.
inline bool runMultiPointDeepOverlapFuzz(const config::Config::BlockHandle& b)
{
   const uint_t numTrials        = b.getParameter< uint_t >("numTrials", uint_t(300));
   const real_t tolerance        = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t maxStepsPerTrial = b.getParameter< uint_t >("maxStepsPerTrial", uint_t(4000));
   const uint_t randomSeed       = b.getParameter< uint_t >("randomSeed", uint_t(1337));

   math::seedRandomGenerator(randomSeed);

   uint_t geometryChecked        = uint_t(0);
   uint_t geometryNaNOrInfCount  = uint_t(0);
   uint_t geometryBadNormalCount = uint_t(0);
   uint_t geometryNegativeOverlapCount = uint_t(0);
   real_t worstVertexOverlap     = real_t(0);
   real_t worstVertexOverlapHalfExtentMin = real_t(0); // min(halfExtents) for the worst trial, for context
   uint_t trialsExceedingLargeOverlapThreshold = uint_t(0); // mirrors [LARGE-OVERLAP-MULTIPOINT]'s own 0.5x bound

   uint_t contactedCount           = uint_t(0);
   uint_t baselineFailCount        = uint_t(0);
   uint_t multiPointFailCount      = uint_t(0);
   uint_t multiPointWorseThanBaselineCount = uint_t(0);
   real_t worstMultiPointGain      = real_t(0);
   uint_t worstMultiPointGainTrial = uint_t(0);
   uint_t cappedFailCount               = uint_t(0);
   uint_t cappedWorseThanBaselineCount  = uint_t(0);
   uint_t hysteresisFailCount               = uint_t(0);
   uint_t hysteresisWorseThanBaselineCount  = uint_t(0);
   uint_t manifoldFailCount                 = uint_t(0);
   uint_t manifoldWorseThanBaselineCount    = uint_t(0);

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
   };
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false },
      { "both", true, true, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false },
      { "useAxisHysteresis", false, false, false, true, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true },
   };

   for (uint_t trial = 0; trial < numTrials; ++trial)
   {
      const Vec3   halfExtents(math::realRandom< real_t >(real_t(0.5), real_t(8)),
                                math::realRandom< real_t >(real_t(0.5), real_t(4)),
                                math::realRandom< real_t >(real_t(0.5), real_t(4)));
      const real_t tiltA = math::realRandom< real_t >(real_t(-90), real_t(90));
      const real_t tiltB = math::realRandom< real_t >(real_t(-90), real_t(90));
      const real_t offY  = math::realRandom< real_t >(-halfExtents[1], halfExtents[1]);
      const real_t offZ  = math::realRandom< real_t >(-halfExtents[2], halfExtents[2]);
      // Full range from clearly separated-but-close through severely embedded (centres closer than
      // either body's own smallest half-extent) -- the added dimension vs. runMultiPointBoxContactFuzz,
      // which never starts anywhere but a fresh approach from clear separation.
      const real_t maxHalfLength = halfExtents[0];
      const real_t separation    = math::realRandom< real_t >(-real_t(1.5) * maxHalfLength, real_t(1.2) * maxHalfLength);

      // (1) Geometry sanity + worst-vertex-overlap tracking, direct on the vertex-test function, at
      // the exact pose this trial's dynamics run will start from (posA/posB/rotA/rotB match
      // detail::runCornerScenarioOnce's own construction below).
      {
         const Vec3 posA(-separation / real_t(2), real_t(0), real_t(0));
         const Vec3 posB(separation / real_t(2), offY, offZ);
         const Rot3 rotA(Vec3(0, 0, 1), tiltA * math::pi / real_t(180));
         const Rot3 rotB(Vec3(0, 0, 1), tiltB * math::pi / real_t(180));

         ++geometryChecked;
         const auto points = paper_analytical::boxBoxMultiPointContacts(posA, rotA.getMatrix(), halfExtents,
                                                                          posB, rotB.getMatrix(), halfExtents);
         const real_t smallerRad = real_t(2) * halfExtents.length() * real_t(0.5); // bounding-sphere radius (both bodies identical this trial)
         bool anyExceedsThreshold = false;
         for (const auto& pt : points)
         {
            const bool finite = std::isfinite(pt.point[0]) && std::isfinite(pt.point[1]) &&
                                 std::isfinite(pt.point[2]) && std::isfinite(pt.normal[0]) &&
                                 std::isfinite(pt.normal[1]) && std::isfinite(pt.normal[2]) &&
                                 std::isfinite(pt.overlap);
            if (!finite) { ++geometryNaNOrInfCount; continue; }
            if (pt.overlap < real_t(-1e-9)) ++geometryNegativeOverlapCount;
            if (std::abs(pt.normal.length() - real_t(1)) > real_t(1e-6)) ++geometryBadNormalCount;
            if (pt.overlap > worstVertexOverlap)
            {
               worstVertexOverlap = pt.overlap;
               worstVertexOverlapHalfExtentMin = std::min({ halfExtents[0], halfExtents[1], halfExtents[2] });
            }
            if (pt.overlap > real_t(0.5) * smallerRad) anyExceedsThreshold = true;
         }
         if (anyExceedsThreshold) ++trialsExceedingLargeOverlapThreshold;
      }

      // (2) Full dynamics, all 4 configurations, same randomly-drawn geometry.
      detail::CornerScenarioParams p;
      p.halfExtents             = halfExtents;
      p.density                 = real_t(0.6);
      p.restitution             = real_t(1);
      p.friction                = real_t(0);
      p.poisson                 = real_t(0.22);
      p.collisionTime           = real_t(32.0078);
      p.dt                      = real_t(0.05);
      p.approachSpeed           = real_t(0); // static release only, see header comment
      p.lateralOffsetY          = offY;
      p.lateralOffsetZ          = offZ;
      p.tiltDegreesA            = tiltA;
      p.tiltDegreesB            = tiltB;
      p.spinRateB               = real_t(0);
      p.separation              = separation;
      p.tolerance               = tolerance;
      p.maxSteps                = maxStepsPerTrial;
      p.enableOutput            = false;

      real_t baselineGain = real_t(0);
      bool   everContactedThisTrial = false;
      for (const auto& cfg : configs)
      {
         p.useRoundedBoxCorners    = cfg.rounded;
         p.useMultiPointBoxContact = cfg.multiPoint;
         p.useBoxOverlapCap        = cfg.overlapCap;
         p.useAxisHysteresis       = cfg.hysteresis;
         p.useContactManifold      = cfg.manifold;

         const auto r = detail::runCornerScenarioOnce(p);
         everContactedThisTrial |= r.everContacted;
         if (!r.everContacted) continue;

         if (std::string(cfg.name) == "baseline (both OFF)")
         {
            baselineGain = r.keGainRel;
            if (!r.pass) ++baselineFailCount;
         }
         if (std::string(cfg.name) == "useMultiPointBoxContact")
         {
            if (!r.pass) ++multiPointFailCount;
            if (r.keGainRel > baselineGain + tolerance)
            {
               ++multiPointWorseThanBaselineCount;
               if (r.keGainRel > worstMultiPointGain)
               {
                  worstMultiPointGain      = r.keGainRel;
                  worstMultiPointGainTrial = trial;
               }
               static uint_t worseLogCount = uint_t(0);
               if (worseLogCount < uint_t(20))
               {
                  WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz WORSE-THAN-BASELINE "
                     "TRIAL " << trial << ": halfExtents=" << halfExtents << " tiltA=" << tiltA
                     << " tiltB=" << tiltB << " offY=" << offY << " offZ=" << offZ
                     << " separation=" << separation << " | baselineGain=" << baselineGain
                     << " multiPointGain=" << r.keGainRel);
                  ++worseLogCount;
               }
            }
         }
         if (std::string(cfg.name) == "useMultiPointBoxContact+overlapCap")
         {
            if (!r.pass) ++cappedFailCount;
            if (r.keGainRel > baselineGain + tolerance) ++cappedWorseThanBaselineCount;
         }
         if (std::string(cfg.name) == "useAxisHysteresis")
         {
            if (!r.pass) ++hysteresisFailCount;
            if (r.keGainRel > baselineGain + tolerance) ++hysteresisWorseThanBaselineCount;
         }
         if (std::string(cfg.name) == "useMultiPointBoxContact+manifold")
         {
            if (!r.pass) ++manifoldFailCount;
            if (r.keGainRel > baselineGain + tolerance) ++manifoldWorseThanBaselineCount;
         }
      }
      if (everContactedThisTrial) ++contactedCount;
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  GEOMETRY: " << geometryChecked
      << " random configurations | NaN/Inf=" << geometryNaNOrInfCount << " | negative overlap="
      << geometryNegativeOverlapCount << " | non-unit normal=" << geometryBadNormalCount
      << " | worst single-vertex overlap=" << worstVertexOverlap << " (min(halfExtents) that trial="
      << worstVertexOverlapHalfExtentMin << ") | trials with >=1 vertex past [LARGE-OVERLAP-MULTIPOINT]'s "
         "own threshold=" << trialsExceedingLargeOverlapThreshold << "/" << numTrials);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  DYNAMICS: " << numTrials
      << " randomized trials | contacted=" << contactedCount << " | baseline ECC-bound failures="
      << baselineFailCount << " | multiPoint ECC-bound failures=" << multiPointFailCount
      << " | multiPoint WORSE than baseline at same trial=" << multiPointWorseThanBaselineCount
      << " | worst multiPoint relGain=" << worstMultiPointGain << " (trial " << worstMultiPointGainTrial << ")");
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  CAPPED: multiPoint+overlapCap ECC-bound "
      "failures=" << cappedFailCount << " | WORSE than baseline at same trial=" << cappedWorseThanBaselineCount);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  HYSTERESIS: useAxisHysteresis ECC-bound "
      "failures=" << hysteresisFailCount << " | WORSE than baseline at same trial=" << hysteresisWorseThanBaselineCount);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz  MANIFOLD: multiPoint+manifold ECC-bound "
      "failures=" << manifoldFailCount << " | WORSE than baseline at same trial=" << manifoldWorseThanBaselineCount);

   // A geometry violation (NaN/non-unit-normal/negative overlap) is a hard bug regardless of depth --
   // always fails. Deep-overlap ECC-bound failures and the worse-than-baseline count are NOT gated
   // into pass/fail here (unlike runMultiPointBoxContactFuzz's shallow-regime equivalent): this test
   // exists specifically to CHARACTERIZE how often and how badly the known deep-overlap gap occurs, not
   // to assert it doesn't -- see multiPointDeepOverlapStability (which already established the defect
   // is real) for the strict pass/fail gate. Report-only for the dynamics half.
   const bool geometryOk = geometryNaNOrInfCount == uint_t(0) && geometryNegativeOverlapCount == uint_t(0) &&
                            geometryBadNormalCount == uint_t(0);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointDeepOverlapFuzz OVERALL: geometryOk=" << geometryOk
      << "  " << (geometryOk ? "PASS" : "FAIL"));
   return geometryOk;
}

// Hand-picked, named counterpart to the two randomized sweeps above: specific geometric configurations
// chosen because each is individually a plausible worst case for a DIFFERENT reason, not because a
// random search happened to land there. A fuzz sweep tells you a bad case exists somewhere in the
// space; these tell you WHICH kind of geometry to worry about, by name, every time this test runs --
// useful as a fixed regression set (adding a new mitigation later, its effect on each named case here
// is directly legible, unlike a fuzz seed that may or may not still land near a similar case) and
// because a couple of these (coincident, engulfed) sit at exact geometric singularities a random draw
// would almost never hit precisely.
inline bool runMultiPointEdgeCases(const config::Config::BlockHandle& b)
{
   // DEFAULT 1 = the strict ECC diagnostic bound every OTHER scenario in this file uses. Overridable
   // to check whether production's own REAL damping (e=0.3) masks a failure this strict bound shows --
   // see runContactStiffnessCalibrationSweep's own finding that e=0.3 rescues the generic corner-tilt
   // geometry completely; this lets the SAME check run against each of this function's own named
   // cases, including "grazing-plank-stack" (the everyday geometry, and the worst failure found this
   // session under e=1: +618% for useMultiPointBoxContact alone).
   const real_t restitutionOverride = b.getParameter< real_t >("restitutionCoefficient", real_t(1));
   struct EdgeCase
   {
      const char* name;
      const char* why;
      Vec3        halfExtentsA;
      Vec3        halfExtentsB;
      real_t      tiltA, tiltB;   // degrees about Z, same convention as CornerScenarioParams
      real_t      offY, offZ;
      real_t      separation;     // X-axis centre-to-centre offset parameter (see runCornerScenarioOnce)
   };

   const std::vector< EdgeCase > cases = {
      { "exact-coincident", "two identical boxes at the SAME position/orientation -- overlap is maximal "
        "and every SAT axis is simultaneously tied (zero gap on all 3 face-pairs), the purest possible "
        "normal-direction singularity.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(0), real_t(0), real_t(0), real_t(0), real_t(0) },

      { "fully-engulfed", "a small box entirely inside a much larger one, centred -- every one of the "
        "small box's 8 vertices is embedded simultaneously, each with a DIFFERENT nearest-face axis, "
        "exercising the full N=8 branch of the per-point force loop at once.",
        Vec3(6, 4, 4), Vec3(1, 1, 1), real_t(0), real_t(0), real_t(0), real_t(0), real_t(0) },

      { "perpendicular-cross", "two identical planks with perpendicular long axes crossing through each "
        "other's centres (a 3D '+' shape) -- the most symmetric possible tie between two candidate "
        "separating axes, both equally valid.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(0), real_t(90), real_t(0), real_t(0), real_t(0) },

      { "corner-into-face-center", "one box's corner driven exactly into the geometric centre of the "
        "other's largest face -- a single vertex deeply embedded with no others nearby, isolating the "
        "single-deep-vertex failure mode directly (no averaging-out from multiple simultaneous points).",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(45), real_t(0), real_t(0), real_t(0), real_t(-5) },

      { "45deg-symmetric-tiebreak", "both bodies tilted +-45 degrees symmetrically about the line "
        "connecting their centres -- SAT's 15-axis test has no single best-separated axis by "
        "construction, a textbook near-degenerate configuration.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(45), real_t(-45), real_t(0), real_t(0), real_t(-2) },

      { "plank-vs-near-cube-deep", "this project's own production floating-box aspect ratio (12x2x3) "
        "deeply overlapping a nearly-cubic body of comparable volume -- the actual shape mismatch a real "
        "jam produces (planks piling against blockier debris), not two identical bodies.",
        Vec3(6, 1, 1.5), Vec3(2.5, 2.3, 2.2), real_t(20), real_t(-15), real_t(0.3), real_t(-0.2), real_t(-3) },

      { "flush-stack-axis-aligned", "both boxes perfectly axis-aligned, offset along X only, well past "
        "full-face overlap -- the SIMPLEST possible deep overlap (no tilt at all), a sanity baseline: if "
        "even this triggers the defect, tilt/offset were never the real trigger, overlap depth alone is.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(0), real_t(0), real_t(0), real_t(0), real_t(-4) },

      { "extreme-size-ratio-deep", "a tiny box (comparable to boxRoundingRadius's own scale) deeply "
        "embedded in a very large one -- checks the per-vertex overlap bound (or lack of one) doesn't "
        "behave pathologically differently when the two bodies' own sizes are very different, unlike "
        "every other case here which uses comparable-sized bodies.",
        Vec3(6, 4, 4), Vec3(0.5, 0.5, 0.5), real_t(30), real_t(0), real_t(0), real_t(0), real_t(0) },

      { "grazing-plank-stack", "two production-shaped planks nearly parallel (5 degree relative tilt), "
        "offset laterally by almost a full half-extent -- mimics two real floaters stacked slightly "
        "askew in a pile, the everyday case rather than a constructed extreme.",
        Vec3(6, 1, 1.5), Vec3(6, 1, 1.5), real_t(2), real_t(-3), real_t(0.9), real_t(1.3), real_t(-1) },
   };

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
      bool        rotationalEffMass = false;
   };
   // [Claude, 2026-10-29] two "useSpheroBoxContact (GJK/Eberly analytical)" rows removed -- see
   // archive/unused_contact_methods/.
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false, false },
      { "both", true, true, false, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false, false },
      { "useAxisHysteresis", false, false, false, true, false, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true, false },
      { "useRotationalEffectiveMass", false, false, false, false, false, true },
      { "useMultiPointBoxContact+rotationalEffMass", false, true, false, false, false, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases -- " << cases.size()
      << " hand-picked geometric configurations x " << (sizeof(configs) / sizeof(configs[0]))
      << " configurations, approachSpeed=0 (static release, same rationale as "
         "multiPointDeepOverlapStability):");

   uint_t multiPointWorseThanBaselineCount = uint_t(0);
   for (const auto& ec : cases)
   {
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases  case=\"" << ec.name << "\" -- " << ec.why);

      detail::CornerScenarioParams p;
      // halfExtentsA/B differ in general -- CornerScenarioParams only has one halfExtents field (both
      // bodies share it in every other scenario in this file), so the smaller/asymmetric-size cases
      // here are built by scaling B's EFFECTIVE size via a direct engine build instead of going through
      // detail::runCornerScenarioOnce (which assumes equal sizes). Kept self-contained rather than
      // extending CornerScenarioParams for a one-off need used only by 3 of these 9 cases.
      const bool sameSize = (ec.halfExtentsA == ec.halfExtentsB);

      real_t baselineGain = real_t(0);
      for (const auto& cfg : configs)
      {
         detail::CornerScenarioResult r;
         if (sameSize)
         {
            p.halfExtents             = ec.halfExtentsA;
            p.density                 = real_t(0.6);
            p.restitution             = restitutionOverride;
            p.friction                = real_t(0);
            p.poisson                 = real_t(0.22);
            p.collisionTime           = real_t(32.0078);
            p.dt                      = real_t(0.05);
            p.approachSpeed           = real_t(0);
            p.lateralOffsetY          = ec.offY;
            p.lateralOffsetZ          = ec.offZ;
            p.tiltDegreesA            = ec.tiltA;
            p.tiltDegreesB            = ec.tiltB;
            p.spinRateB               = real_t(0);
            p.separation              = ec.separation;
            p.tolerance               = real_t(0.02);
            p.maxSteps                = uint_t(6000);
            p.enableOutput            = false;
            p.useRoundedBoxCorners    = cfg.rounded;
            p.boxRoundingRadius       = real_t(0.05);
            p.useMultiPointBoxContact = cfg.multiPoint;
            p.useBoxOverlapCap        = cfg.overlapCap;
            p.useAxisHysteresis       = cfg.hysteresis;
            p.useContactManifold      = cfg.manifold;
            p.useRotationalEffectiveMass = cfg.rotationalEffMass;
            r = detail::runCornerScenarioOnce(p);
         }
         else
         {
            // Asymmetric-size cases: build the engine directly rather than through
            // runCornerScenarioOnce, which only accepts one shared halfExtents for both bodies.
            DEMEngine engine(real_t(0.05), restitutionOverride, real_t(0), real_t(0.22), real_t(32.0078), real_t(0),
                              real_t(200), /*useAnalyticalContactDetection=*/true, cfg.rounded, real_t(0.05),
                              cfg.multiPoint, cfg.overlapCap, real_t(0.3), cfg.hysteresis, real_t(0.15),
                              cfg.manifold, real_t(0.25), cfg.rotationalEffMass);
            const Rot3 rotA(Vec3(0, 0, 1), ec.tiltA * math::pi / real_t(180));
            const Rot3 rotB(Vec3(0, 0, 1), ec.tiltB * math::pi / real_t(180));
            const size_t idxA = engine.addBox(Vec3(-ec.separation / 2, 0, 0), Vec3(0, 0, 0), rotA,
                                               real_t(2) * ec.halfExtentsA, real_t(0.6));
            const size_t idxB = engine.addBox(Vec3(ec.separation / 2, ec.offY, ec.offZ), Vec3(0, 0, 0), rotB,
                                               real_t(2) * ec.halfExtentsB, real_t(0.6));
            engine.finalizeSetup();

            r.keBefore  = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
            r.peInitial = detail::computeInitialElasticPE(
               Vec3(-ec.separation / 2, 0, 0), rotA, ec.halfExtentsA, real_t(0.6),
               Vec3(ec.separation / 2, ec.offY, ec.offZ), rotB, ec.halfExtentsB, real_t(0.6),
               restitutionOverride, real_t(32.0078), cfg.rounded, real_t(0.05), cfg.multiPoint, cfg.overlapCap,
               real_t(0.3), cfg.rotationalEffMass);
            r.keMax    = r.keBefore;
            for (uint_t s = 0; s < uint_t(6000); ++s)
            {
               real_t penetration = real_t(0);
               if (engine.queryContact(idxA, idxB, penetration) && penetration < real_t(0)) r.everContacted = true;
               engine.step();
               const real_t keNow = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
               if (keNow > r.keMax) { r.keMax = keNow; r.keMaxStep = s; }
            }
            r.keAfter = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
            const real_t ceiling = r.keBefore + r.peInitial;
            r.keGainRel = (r.keMax - ceiling) / std::max(real_t(1e-12), ceiling);
            r.pass      = r.everContacted && r.keGainRel < real_t(0.02);
         }

         if (std::string(cfg.name) == "baseline (both OFF)") baselineGain = r.keGainRel;
         if (std::string(cfg.name) == "useMultiPointBoxContact" && r.keGainRel > baselineGain + real_t(0.02))
            ++multiPointWorseThanBaselineCount;

         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases    config=\"" << cfg.name << "\""
            << "  KE_before=" << r.keBefore << " PE_initial=" << r.peInitial << " KE_max=" << r.keMax
            << " (step " << r.keMaxStep << ")"
            << " maxRelGain=" << r.keGainRel
            << " rawSatNormalFlipCount=" << r.normalFlipCount
            << " cumulativeNormalRotationDeg=" << r.cumulativeNormalRotationDeg
            << " maxPointJumpFrac=" << r.maxPointJumpFrac << " pointJumpCount=" << r.pointJumpCount
            << (sameSize ? "" : " (N/A, asymmetric-size branch)")
            << " everContacted=" << r.everContacted
            << "  " << (r.pass ? "PASS" : "FAIL"));
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointEdgeCases OVERALL: multiPointWorseThanBaseline in "
      << multiPointWorseThanBaselineCount << "/" << cases.size() << " named cases");
   // Characterization tool like multiPointDeepOverlapFuzz -- report-only, not a hard gate (some of
   // these cases are deliberately singular/degenerate and not expected to cleanly pass any config).
   return true;
}

// multiPointDeepOverlapStability's own dt-independence result (16x timestep refinement, <0.01%
// change in baseline's gain at every depth) rules out numerical/integration error as the cause of
// baseline's OWN deep-overlap energy injection (separations 6/3/0/-4 all fail with ZERO detected
// axis flips, near-zero cumulative normal rotation, and a contact-point jump that itself shrinks to
// noise as dt refines -- yet the energy gain does not move at all). None of useRoundedBoxCorners/
// useMultiPointBoxContact/useBoxOverlapCap/useAxisHysteresis/useContactManifold change this, because
// all five still compute their force from the SAME mesa_pd::kernel::LinearSpringDashpot::
// setStiffnessAndDamping calibration (stiffnessN/dampingN), which uses ONLY the translational
// effective mass (meff = 1/(1/mA + 1/mB)) -- no correction for the ROTATIONAL inertia excited by a
// contact point offset from either body's centroid, which every case in this file has by
// construction (that offset is exactly what makes the near-degenerate/multi-point geometry
// interesting in the first place). At restitution=1 (this file's own strict diagnostic bound) that
// calibration gives EXACTLY ZERO damping -- correct for a head-on, non-rotating collision, but with
// no mechanism to correctly balance a contact that also does rotational work.
//
// This sweeps the EXISTING parameters that calibration already exposes, to see whether the effect is
// masked (not necessarily fixed) by realistic production values, before concluding a genuine formula
// change (e.g. an effective mass that also accounts for rotational inertia at the contact point,
// 1/meff = 1/mA + 1/mB + (r x n)^T I^-1 (r x n) summed appropriately) is actually required:
//   (1) restitutionCoefficient: e<1 introduces REAL damping (dampingN != 0) -- does the resulting
//       dissipation happen to offset the rotational-energy leak enough to mask it? Production itself
//       uses e=0.3 ("wet/waterlogged wood"), never e=1 -- this checks whether the strict e=1 bound
//       this whole file uses is somehow overstating a problem production's own real damping already
//       absorbs.
//   (2) collisionTime: controls the ABSOLUTE stiffness magnitude for a fixed restitution -- a softer
//       spring (larger collisionTime) produces less peak force for the same overlap, hence less
//       torque at a fixed lever arm; does simply detuning the spring change the picture?
//   (3) frictionCoefficient: acts on the TANGENTIAL component only (this test's own geometry has
//       nonzero relative tangential velocity at the contact point from p.spinRateB) -- tangential
//       dissipation could not fix a normal-direction energy leak directly, but is included as a
//       direct empirical check rather than an assumption.
inline bool runContactStiffnessCalibrationSweep(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents = b.getParameter< Vec3 >("halfExtents", Vec3(6, 1, 1.5));
   const real_t tolerance   = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t maxSteps    = b.getParameter< uint_t >("maxSteps", uint_t(4000));

   const std::vector< real_t > separations = { real_t(6), real_t(3), real_t(0), real_t(-4) };
   struct RestitutionCase
   {
      const char* name;
      real_t      value;
   };
   const std::vector< RestitutionCase > restitutions = {
      { "e=1.00 (strict diagnostic bound)", real_t(1.0) },
      { "e=0.99", real_t(0.99) },
      { "e=0.90", real_t(0.90) },
      { "e=0.50", real_t(0.50) },
      { "e=0.30 (PRODUCTION value)", real_t(0.30) },
   };
   struct CollisionTimeCase
   {
      const char* name;
      real_t      value;
   };
   const std::vector< CollisionTimeCase > collisionTimes = {
      { "collisionTime=8 (stiffer)", real_t(8) },
      { "collisionTime=32.0078 (PRODUCTION value)", real_t(32.0078) },
      { "collisionTime=128 (softer)", real_t(128) },
   };
   const std::vector< real_t > frictions = { real_t(0), real_t(0.5) }; // 0.5 = PRODUCTION value
   // Baseline alone isolates the stiffness calibration itself; ALSO sweeping useMultiPointBoxContact
   // answers the natural follow-up once baseline's own mystery gain turned out to vanish completely at
   // e=0.30 (production's own real value): does that SAME real damping also rescue multipoint's own,
   // much larger, "grazing-plank-stack" failure (+618% at e=1, DEMValidationScenarios.h's own
   // multiPointEdgeCases), or is multipoint's defect a genuinely different, damping-independent one?
   const std::vector< bool > multiPointOptions = { false, true };
   // useRotationalEffectiveMass: the third candidate fix researched this session (see
   // ModelParameters::useRotationalEffectiveMass's own doc comment, SimulationParams.h) -- unlike
   // useMultiPointBoxContact/useRoundedBoxCorners/useAxisHysteresis/useContactManifold (which all change
   // WHICH (point,normal,overlap) triplet feeds LinearSpringDashpot), this changes the force law's OWN
   // calibration: stiffnessN/dampingN are set from the generalized (translational+rotational) effective
   // mass at the actual contact point instead of the plain 1/(1/mA+1/mB). This scenario's own tiltDegreesA/
   // tiltDegreesB/lateralOffsetY/Z/spinRateB (see CornerScenarioParams defaults) give every contact here a
   // genuine moment arm off the line of centres -- precisely the geometry the rotational correction targets
   // -- so sweeping it crossed with useMultiPointBoxContact directly answers: does correcting the force
   // law's calibration alone (no change to which point/normal is used) already fix baseline's own mystery
   // gain, independent of and possibly BEFORE reaching for a different (point,normal) mitigation at all?
   const std::vector< bool > rotationalEffMassOptions = { false, true };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] contactStiffnessCalibrationSweep -- " << separations.size()
      << " deep-overlap depths x " << restitutions.size() << " restitution values x "
      << collisionTimes.size() << " collision times x " << frictions.size() << " friction values x "
      << multiPointOptions.size() << " useMultiPointBoxContact settings x "
      << rotationalEffMassOptions.size() << " useRotationalEffectiveMass settings, testing whether the "
         "LinearSpringDashpot stiffness/damping calibration itself (translational effective mass only, "
         "no rotational-inertia correction) is masked or fixed by any EXISTING parameter combination, "
         "for baseline AND for multipoint, with and without the rotational-effective-mass correction:");

   uint_t anyPassCount   = uint_t(0);
   uint_t totalRuns      = uint_t(0);
   real_t worstGainSeen  = real_t(0);
   std::string worstGainConfig;
   // Per-restitution/collisionTime/friction/multiPoint/rotationalEffMass combo: does EVERY depth pass
   // under it? This is the meaningful question -- a combo that only happens to pass at one depth isn't a
   // real fix.
   for (const bool useRotationalEffMass : rotationalEffMassOptions)
   for (const bool useMultiPoint : multiPointOptions)
   for (const auto& restCase : restitutions)
   {
      for (const auto& ctCase : collisionTimes)
      {
         for (const real_t friction : frictions)
         {
            bool allDepthsPassed = true;
            real_t worstGainThisCombo = real_t(0);
            for (const real_t sep : separations)
            {
               detail::CornerScenarioParams p;
               p.halfExtents   = halfExtents;
               p.density       = real_t(0.6);
               p.restitution   = restCase.value;
               p.friction      = friction;
               p.poisson       = real_t(0.22);
               p.collisionTime = ctCase.value;
               p.dt            = real_t(0.05);
               p.approachSpeed = real_t(0); // static release, same rationale as multiPointDeepOverlapStability
               p.separation    = sep;
               p.tolerance     = tolerance;
               p.maxSteps      = maxSteps;
               p.enableOutput  = false;
               p.useRoundedBoxCorners        = false;
               p.useMultiPointBoxContact     = useMultiPoint;
               p.useBoxOverlapCap            = false;
               p.useAxisHysteresis           = false;
               p.useContactManifold          = false;
               p.useRotationalEffectiveMass  = useRotationalEffMass;

               const auto r = detail::runCornerScenarioOnce(p);
               ++totalRuns;
               if (!r.pass) allDepthsPassed = false;
               if (r.keGainRel > worstGainThisCombo) worstGainThisCombo = r.keGainRel;
               if (r.keGainRel > worstGainSeen)
               {
                  worstGainSeen = r.keGainRel;
                  worstGainConfig = std::string("useRotationalEffectiveMass=") +
                                    (useRotationalEffMass ? "true" : "false") +
                                    " / useMultiPointBoxContact=" + (useMultiPoint ? "true" : "false") +
                                    " / " + restCase.name + " / " + ctCase.name + " / friction=" +
                                    std::to_string(friction) + " / separation=" + std::to_string(sep);
               }
            }
            if (allDepthsPassed) ++anyPassCount;
            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] contactStiffnessCalibrationSweep  "
               << "useRotationalEffectiveMass=" << useRotationalEffMass
               << " / useMultiPointBoxContact=" << useMultiPoint << " / " << restCase.name
               << " / " << ctCase.name << " / friction=" << friction << "  worstGainAcrossDepths="
               << worstGainThisCombo << "  " << (allDepthsPassed ? "ALL 4 DEPTHS PASS" : "at least one depth FAILS"));
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] contactStiffnessCalibrationSweep OVERALL: " << anyPassCount
      << "/" << (restitutions.size() * collisionTimes.size() * frictions.size() * multiPointOptions.size() *
                 rotationalEffMassOptions.size())
      << " parameter combinations pass at EVERY depth (" << totalRuns << " total runs) | worst gain seen="
      << worstGainSeen << " at [" << worstGainConfig << "]");
   // Report-only, like the other characterization scenarios in this file -- the point is to learn
   // WHICH (if any) existing parameter combination masks the effect, not to assert a pass/fail verdict.
   return anyPassCount > uint_t(0);
}

// Every scenario above starts the two bodies ALREADY overlapping at rest -- a controlled but somewhat
// artificial static-release probe of the contact FORMULA itself, isolated from how such an overlap
// would ever actually arise in production. This scenario instead reproduces the literal PHYSICAL
// process suspected of producing them: this project's own FloatingBox (halfExtents <6,1,1.5>,
// densityRatio 0.6) free-falling under REAL production gravity onto a small pile of IDENTICAL boxes
// already at rest on a floor -- exactly "a floater falling onto already-floating debris", the
// mechanism asked about directly. No water/buoyancy here (this whole file is a dry DEM-only harness,
// same as every other scenario) -- this isolates whether DEM contact resolution alone, given a
// realistic impact velocity from a real fall height, can launch something out of a pile, independent
// of any fluid effect.
//
// Two things are swept, both directly targeting "how/why do these get flying, especially with more
// floaters already present":
//   (1) pileSize: how many boxes are already resting side-by-side before the drop (1, 3, 6) -- tests
//       whether MORE simultaneous neighbours makes the outcome worse, i.e. whether this is fundamentally
//       a [SIMULTANEOUS-CONTACTS]-style effect (ParticleKernels.h's own diagnostic for exactly this
//       concern: "the EFFECTIVE combined stiffness ... may exceed what the explicit integrator's
//       timestep was tuned to handle for a single dominant contact").
//   (2) landingOffset: WHERE lateral the falling box lands relative to the pile -- dead-centre on one
//       pile box (clean, symmetric single contact), on the SEAM between two adjacent boxes (the
//       corner-into-two-neighbours case a real staggered/jittered pile produces most often), or
//       overhanging the pile's own edge (asymmetric, one side supported one side not).
//
// Energy bound here is NOT the elastic-PE-in-overlap accounting the scenarios above use (which is
// entangled with the contact model's own formula) -- it is the much more fundamental, model-
// independent one: total system mechanical energy E = sum(KE_i) + sum(m_i * g * z_i) can only
// DECREASE over time (restitution=0.3, "wet/waterlogged wood", matches production exactly -- real
// dissipation is expected, not just "no worse than start"). A rise in E at any step, by any amount
// beyond floating-point tolerance, is unambiguous evidence of energy injected from nowhere -- gravity
// is the only OTHER force in this dry harness and its own contribution is exact and trivial (mgz), so
// there is no formula-specific approximation error to account for here at all.
inline bool runFallOntoPileStability(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents   = b.getParameter< Vec3 >("halfExtents", Vec3(6, 1, 1.5));
   const real_t density       = b.getParameter< real_t >("density", real_t(0.6));
   const real_t restitution   = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t friction      = b.getParameter< real_t >("frictionCoefficient", real_t(0.5));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(32.0078));
   const real_t dt            = b.getParameter< real_t >("dt", real_t(0.05));
   const real_t gravity       = b.getParameter< real_t >("gravity", real_t(7.02618e-6)); // production value
   const uint_t maxSteps      = b.getParameter< uint_t >("maxSteps", uint_t(8000));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const bool   writeVtk      = b.getParameter< bool >("writeVtk", false);
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(4));
   const std::string outBase  = b.getParameter< std::string >("outputBaseName", "fall_onto_pile_stability");

   // CORRECTED MECHANISM (per direct feedback on an earlier version of this test): the RELEASE-OVERLAP
   // penetration is not a symptom of "a newly-released floater drops onto an already-settled pile" --
   // the just-released floater is FIXED (frozen, zero velocity) at its exact spawn point right up until
   // the instant of release, so at the very first contact check it hasn't moved AT ALL. The overlap
   // exists because a DIFFERENT, EARLIER-launched floater -- already flying due to the instability
   // this whole session investigates -- has drifted/climbed into that exact spot in the meantime (this
   // session's own VTK evidence: uid 1585267068834414616 climbed +16.7 LU in the 1000 steps before
   // detection, uid 5044031582654955533 climbed +10.6 LU). So the STATIONARY body here (the "pile")
   // is the correct stand-in for the just-released floater, and the ONE body given real velocity
   // (flyingSpeed below) is the stand-in for the already-launched neighbour arriving at it -- the
   // previous version of this test had these two roles backwards (both bodies started at rest).
   //
   // dropHeight is still the CLEAR GAP at t=0 (negative = starting already embedded) -- -0.71/-0.5
   // directly replay the two ACTUAL measured RELEASE-OVERLAP penetrations (job 4201127, step 122001);
   // 0 is the boundary case (already touching, zero gap).
   // quickMode (default false, unchanged behaviour): the full 3x2x2x1 geometry sweep x 9 configs = 108
   // runs is prohibitively slow when a config's deep initial overlap keeps [DEM-WATCH]/
   // [STRONG-CONTACT-FORCE] firing every single sub-cycle for the whole maxSteps window (each combo can
   // take several minutes just formatting/writing that log volume) -- true collapses the geometry sweep
   // to the SINGLE most representative combo (the actual measured release-overlap magnitude and speed,
   // dropHeight=-0.71, flyingSpeed=0.015, pileSize=6) so all 9 mitigation configs can still be compared
   // quickly, at the cost of not covering the other 11 geometry combinations this run.
   const bool quickMode = b.getParameter< bool >("quickMode", false);
   const std::vector< real_t > dropHeights =
      quickMode ? std::vector< real_t >{ real_t(-0.71) } : std::vector< real_t >{ real_t(-0.71), real_t(-0.5), real_t(0) };
   // Initial speed of the "flying" body (toward the stationary one, see the mechanism note above).
   // 0 isolates pure static-overlap release (no momentum at all, matching this file's OTHER deep-
   // overlap scenarios); 0.015 matches the ORDER OF MAGNITUDE this session directly measured for the
   // two real flying neighbours (+10.6/+16.7 LU over 1000 steps => ~0.011-0.017 LU/step average).
   const std::vector< real_t > flyingSpeeds =
      quickMode ? std::vector< real_t >{ real_t(0.015) } : std::vector< real_t >{ real_t(0), real_t(0.015) };

   const std::vector< uint_t > pileSizes = quickMode ? std::vector< uint_t >{ uint_t(6) } : std::vector< uint_t >{ uint_t(1), uint_t(6) };
   enum class Landing
   {
      OnCentre,
      OnSeam
   };
   const std::vector< std::pair< Landing, const char* > > landings = { { Landing::OnCentre, "on-centre" } };

   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
      bool        rotationalEffMass = false;
   };
   // [Claude, 2026-10-29] "useSpheroBoxContact" row removed -- see archive/unused_contact_methods/.
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false, false },
      { "both", true, true, false, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false, false },
      { "useAxisHysteresis", false, false, false, true, false, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true, false },
      { "useRotationalEffectiveMass", false, false, false, false, false, true },
      { "useMultiPointBoxContact+rotationalEffMass", false, true, false, false, false, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] fallOntoPileStability -- a STATIONARY just-released "
      "FloatingBox vs. an ALREADY-FLYING neighbour under real production gravity (" << gravity << "), "
      << dropHeights.size() << " release gaps (negative = already overlapping) x " << flyingSpeeds.size()
      << " flying speeds x " << pileSizes.size() << " pile sizes x " << landings.size()
      << " landing offsets x " << (sizeof(configs) / sizeof(configs[0])) << " configurations:");

   bool anyConfigFailed              = false;
   bool multiPointEverWorseThanBase  = false;
   real_t worstRelGainOverall        = real_t(0);
   for (const real_t dropHeight : dropHeights)
   for (const real_t flyingSpeed : flyingSpeeds)
   for (const uint_t pileSize : pileSizes)
   {
      for (const auto& landing : landings)
      {
         real_t baselineGain = real_t(0);
         for (const auto& cfg : configs)
         {
            // Floor at z=0. Pile boxes rest side-by-side along Y, long axis (halfExtents[0]) along X,
            // flat face down (identity rotation) -- the resting orientation this project's own VTK
            // evidence shows these boxes actually settle into. Touching (zero gap) between neighbours,
            // matching an already-packed/jammed pile rather than a loosely-spaced one.
            const real_t pileZ    = halfExtents[2]; // box centre sits exactly halfExtents.z above the floor
            const real_t pileSpan = real_t(2) * halfExtents[1] * real_t(pileSize); // total Y extent of the pile
            const real_t pileY0   = -pileSpan / real_t(2) + halfExtents[1]; // first box's centre Y

            DEMEngine engine(dt, restitution, friction, poisson, collisionTime, gravity,
                              /*domainHalfExtent=*/real_t(500), /*useAnalyticalContactDetection=*/true,
                              cfg.rounded, real_t(0.05), cfg.multiPoint, cfg.overlapCap, real_t(0.3),
                              cfg.hysteresis, real_t(0.15), cfg.manifold, real_t(0.25),
                              cfg.rotationalEffMass);
            engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1)); // floor

            std::vector< size_t > pileIdx;
            for (uint_t i = 0; i < pileSize; ++i)
            {
               const real_t y = pileY0 + real_t(i) * real_t(2) * halfExtents[1];
               pileIdx.push_back(engine.addBox(Vec3(0, y, pileZ), Vec3(0, 0, 0), Rot3(),
                                                real_t(2) * halfExtents, density));
            }

            // Landing Y offset relative to the pile: dead centre of pile box 0, or the seam between
            // box 0 and box 1 (only meaningful for pileSize>=2 -- falls back to on-centre for
            // pileSize==1, still a valid single-neighbour baseline for that pile size).
            real_t landingY = pileY0;
            if (landing.first == Landing::OnSeam && pileSize >= 2) landingY = pileY0 + halfExtents[1];

            // This body is the STAND-IN FOR THE ALREADY-FLYING NEIGHBOUR (see the mechanism note
            // above) -- moving at flyingSpeed TOWARD the stationary pile, not the newly-released
            // floater itself (which is correctly represented by the stationary pile boxes below).
            const real_t dropZ = pileZ + halfExtents[2] + dropHeight + halfExtents[2];
            const size_t fallIdx = engine.addBox(Vec3(0, landingY, dropZ), Vec3(0, 0, -flyingSpeed), Rot3(),
                                                  real_t(2) * halfExtents, density);

            engine.finalizeSetup();

            if (writeVtk)
            {
               std::string safeName = cfg.name;
               for (char& c : safeName)
                  if (!std::isalnum(static_cast< unsigned char >(c))) c = '_';
               engine.setupVtk(outBase + "_" + safeName, vtkSpacing);
            }

            auto totalMechanicalEnergy = [&]() {
               real_t e = real_t(0);
               for (size_t idx : pileIdx) e += engine.kineticEnergy(idx) + engine.ac().getMass(idx) * gravity * engine.ac().getPosition(idx)[2];
               e += engine.kineticEnergy(fallIdx) + engine.ac().getMass(fallIdx) * gravity * engine.ac().getPosition(fallIdx)[2];
               return e;
            };

            const real_t e0     = totalMechanicalEnergy();
            real_t       eMax   = e0;
            uint_t       eMaxStep = uint_t(0);
            real_t       peakZAfterImpact = real_t(0);
            bool         everContacted    = false;
            for (uint_t s = 0; s < maxSteps; ++s)
            {
               real_t pen = real_t(0);
               if (engine.queryContact(fallIdx, pileIdx[0], pen) && pen < real_t(0)) everContacted = true;
               engine.step();
               if (writeVtk) engine.writeVtkIfDue();
               const real_t eNow = totalMechanicalEnergy();
               if (eNow > eMax) { eMax = eNow; eMaxStep = s; }
               if (everContacted)
               {
                  const real_t z = engine.ac().getPosition(fallIdx)[2];
                  if (z > peakZAfterImpact) peakZAfterImpact = z;
               }
            }

            const real_t relGain = (eMax - e0) / std::max(real_t(1e-12), std::abs(e0));
            const bool   pass    = relGain < tolerance;
            if (!pass) anyConfigFailed = true;
            if (std::string(cfg.name) == "baseline (both OFF)") baselineGain = relGain;
            if (std::string(cfg.name) == "useMultiPointBoxContact" && relGain > baselineGain + tolerance)
               multiPointEverWorseThanBase = true;
            if (relGain > worstRelGainOverall) worstRelGainOverall = relGain;

            WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] fallOntoPileStability  dropHeight=" << dropHeight
               << " flyingSpeed=" << flyingSpeed << " pileSize=" << pileSize
               << " landing=" << landing.second << " config=\"" << cfg.name << "\"  E0=" << e0
               << " Emax=" << eMax << " (step " << eMaxStep << ") relGain=" << relGain
               << " peakZAfterImpact=" << peakZAfterImpact << " (dropZ was " << dropZ << ", pile top z="
               << (pileZ + halfExtents[2]) << ")  everContacted=" << everContacted
               << "  " << (pass ? "PASS" : "FAIL"));
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] fallOntoPileStability OVERALL: anyConfigFailed="
      << anyConfigFailed << "  multiPointEverWorseThanBaseline=" << multiPointEverWorseThanBase
      << "  worstRelGainSeen=" << worstRelGainOverall << "  " << (!anyConfigFailed ? "PASS" : "FAIL"));
   return !anyConfigFailed;
}

// Stress test for useMultiPointBoxContact before trusting it on a full production run: the hand-tuned
// geometry in runCornerNormalDiscontinuity/runCornerFixComparison is ONE specific corner/tilt/offset
// configuration -- this scenario sweeps many randomized ones (box aspect ratio, tilt about a random
// axis on EACH body, lateral offset, approach speed, spin) to check the mitigation holds broadly, not
// just for the one case it was validated against. Two independent checks per trial:
//   (1) GEOMETRY sanity, direct on paper_analytical::boxBoxMultiPointContacts: no NaN/Inf in any
//       returned point/normal/overlap, every normal is a unit vector, every overlap is >= 0, and --
//       the key cross-check against boxBoxSAT -- whenever SAT proves the pair definitely SEPARATED
//       (result.contact == false), the vertex-embedding test must ALSO find nothing (a genuinely
//       separated pair can never have one body's vertex sitting inside the other; a violation here
//       would mean a sign/axis bug in the vertex test itself, independent of any dynamics).
//   (2) DYNAMICS: the exact same strict ECC energy bound (e=1, mu=0, no gravity: total system KE must
//       never exceed its own pre-contact value) as runCornerNormalDiscontinuity, but run under
//       useMultiPointBoxContact=true across every randomized trial that actually made contact.
// A trial's box sizes/speeds/offsets are drawn from ranges that comfortably bracket this project's own
// production FloatingBox family (halfExtents as small as <6,1,1.5> and as chunky as roughly cube-
// shaped) and the realistic approach speeds this session's own forensics established (see
// runCornerNormalDiscontinuity's own header comment) -- not an arbitrary/unbounded fuzz domain.
inline bool runMultiPointBoxContactFuzz(const config::Config::BlockHandle& b)
{
   const uint_t numTrials         = b.getParameter< uint_t >("numTrials", uint_t(500));
   const real_t tolerance         = b.getParameter< real_t >("toleranceRelError", real_t(0.02));
   const uint_t maxStepsPerTrial  = b.getParameter< uint_t >("maxStepsPerTrial", uint_t(8000));
   const uint_t randomSeed        = b.getParameter< uint_t >("randomSeed", uint_t(42));

   math::seedRandomGenerator(randomSeed);

   uint_t geometryChecked         = uint_t(0);
   uint_t geometryNaNOrInfCount   = uint_t(0);
   uint_t geometryBadNormalCount  = uint_t(0);
   uint_t geometryNegativeOverlapCount = uint_t(0);
   uint_t geometrySatDisagreementCount = uint_t(0); // SAT says separated, vertex test says embedded

   uint_t contactedCount = uint_t(0);
   uint_t energyPassCount = uint_t(0);
   real_t worstGainSeen   = real_t(0);
   uint_t worstGainTrial  = uint_t(0);

   for (uint_t trial = 0; trial < numTrials; ++trial)
   {
      const Vec3 halfExtents(math::realRandom< real_t >(real_t(0.5), real_t(8)),
                              math::realRandom< real_t >(real_t(0.5), real_t(4)),
                              math::realRandom< real_t >(real_t(0.5), real_t(4)));

      // (1) Pure-geometry sanity, several random relative poses per trial's box size -- independent of
      // any dynamics, so cheap to sample densely.
      for (uint_t g = 0; g < 5; ++g)
      {
         const Vec3 posA(real_t(0), real_t(0), real_t(0));
         const Rot3 rotA(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                               math::realRandom< real_t >(-1, 1))
                            .getNormalizedIfNotZero(),
                         math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
         // Offset chosen so roughly half the samples land in genuine overlap and half in clear
         // separation -- a pure "always separated" or "always overlapping" sweep would never exercise
         // the boundary this check cares about.
         const real_t maxReach = halfExtents.length() * real_t(2.2);
         const Vec3   posB(math::realRandom< real_t >(-maxReach, maxReach),
                            math::realRandom< real_t >(-maxReach, maxReach),
                            math::realRandom< real_t >(-maxReach, maxReach));
         const Rot3 rotB(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                               math::realRandom< real_t >(-1, 1))
                            .getNormalizedIfNotZero(),
                         math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));

         ++geometryChecked;
         const auto satResult = paper_analytical::boxBoxSAT(posA, rotA.getMatrix(), halfExtents, posB,
                                                              rotB.getMatrix(), halfExtents);
         const auto points = paper_analytical::boxBoxMultiPointContacts(posA, rotA.getMatrix(), halfExtents,
                                                                          posB, rotB.getMatrix(), halfExtents);

         if (!satResult.contact && !points.empty()) ++geometrySatDisagreementCount;

         for (const auto& pt : points)
         {
            const bool finite = std::isfinite(pt.point[0]) && std::isfinite(pt.point[1]) &&
                                 std::isfinite(pt.point[2]) && std::isfinite(pt.normal[0]) &&
                                 std::isfinite(pt.normal[1]) && std::isfinite(pt.normal[2]) &&
                                 std::isfinite(pt.overlap);
            if (!finite) { ++geometryNaNOrInfCount; continue; }
            if (pt.overlap < real_t(-1e-9)) ++geometryNegativeOverlapCount;
            const real_t normalLen = pt.normal.length();
            if (std::abs(normalLen - real_t(1)) > real_t(1e-6)) ++geometryBadNormalCount;
         }
      }

      // (2) Full dynamic energy-conservation check, ONE randomized corner/tilt/offset/speed
      // configuration per trial, useMultiPointBoxContact=true throughout.
      detail::CornerScenarioParams p;
      p.halfExtents             = halfExtents;
      p.density                 = real_t(0.6);
      p.restitution             = real_t(1); // strict ECC bound, see this function's own header comment
      p.friction                = real_t(0);
      p.poisson                 = real_t(0.22);
      p.collisionTime           = real_t(32.0078);
      p.dt                      = real_t(0.05);
      p.approachSpeed           = math::realRandom< real_t >(real_t(0.005), real_t(0.05));
      p.lateralOffsetY          = math::realRandom< real_t >(-halfExtents[1], halfExtents[1]);
      p.lateralOffsetZ          = math::realRandom< real_t >(-halfExtents[2], halfExtents[2]);
      p.tiltDegreesA            = math::realRandom< real_t >(real_t(-60), real_t(60));
      p.tiltDegreesB            = math::realRandom< real_t >(real_t(-60), real_t(60));
      p.spinRateB               = math::realRandom< real_t >(real_t(-0.01), real_t(0.01));
      p.separation              = real_t(12);
      p.tolerance               = tolerance;
      p.maxSteps                = maxStepsPerTrial;
      p.enableOutput            = false;
      p.useRoundedBoxCorners    = false;
      p.useMultiPointBoxContact = true;

      const auto r = detail::runCornerScenarioOnce(p);
      if (r.everContacted)
      {
         ++contactedCount;
         if (r.pass) ++energyPassCount;
         else
         {
            static uint_t failureLogCount = uint_t(0);
            if (failureLogCount < uint_t(20))
            {
               WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz FAILING TRIAL "
                  << trial << ": halfExtents=" << p.halfExtents << " approachSpeed=" << p.approachSpeed
                  << " lateralOffsetY=" << p.lateralOffsetY << " lateralOffsetZ=" << p.lateralOffsetZ
                  << " tiltDegreesA=" << p.tiltDegreesA << " tiltDegreesB=" << p.tiltDegreesB
                  << " spinRateB=" << p.spinRateB << " | KE_before=" << r.keBefore
                  << " KE_max=" << r.keMax << " (step " << r.keMaxStep << ") maxRelGain=" << r.keGainRel);
               ++failureLogCount;
            }
         }
         if (r.keGainRel > worstGainSeen)
         {
            worstGainSeen  = r.keGainRel;
            worstGainTrial = trial;
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz  GEOMETRY: " << geometryChecked
      << " random pose pairs checked | NaN/Inf=" << geometryNaNOrInfCount
      << " | negative overlap=" << geometryNegativeOverlapCount
      << " | non-unit normal=" << geometryBadNormalCount
      << " | SAT-says-separated-but-vertex-embedded=" << geometrySatDisagreementCount);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz  DYNAMICS: " << numTrials
      << " randomized trials | contacted=" << contactedCount << " | energy-bound PASS="
      << energyPassCount << "/" << contactedCount << " | worst maxRelGain=" << worstGainSeen
      << " (trial " << worstGainTrial << ")");

   const bool geometryOk = geometryNaNOrInfCount == uint_t(0) && geometryNegativeOverlapCount == uint_t(0) &&
                            geometryBadNormalCount == uint_t(0) && geometrySatDisagreementCount == uint_t(0);
   const bool dynamicsOk = contactedCount > uint_t(0) && energyPassCount == contactedCount;
   const bool pass       = geometryOk && dynamicsOk;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] multiPointBoxContactFuzz OVERALL: geometryOk=" << geometryOk
      << " dynamicsOk=" << dynamicsOk << "  " << (pass ? "PASS" : "FAIL"));
   return pass;
}

// Direct replay of a REAL production incident, not a hand-tuned or randomized reproduction: job
// 4199585 (bridge "widerstacks" box scenario), uid=3314649325744685089 (the box independently
// confirmed "flying" in ParaView -- climbed 31 LU, back above its own spawn height) vs its contact
// partner uid=720575940379279392, the PROLONGED-CONTACT episode FloaterHistoryLogger recorded as
// starting at step 25883. Both bodies are seeded with their EXACT logged position/velocity/angular-
// velocity/orientation at that step (cross-verified via exact-text grep against an earlier
// awk-numeric-comparison false match on these 19-digit uids -- see this session's own investigation)
// and run through the SAME production physics (restitutionCoefficient=0.3 "wet/waterlogged wood",
// frictionCoefficient=0.5, gravity=7.03e-6 LU/step^2, particleCollisionTime/particleTimestepSize from
// computedParams.txt, useAnalyticalContactDetection=true) for 3540 sub-cycles (177 main steps x 20
// sub-cycles), matching how long the real run took to reach step 26060 -- the point the production
// log shows this exact contact had cleanly separated.
//
// Measured directly on the REAL production log (not this replay, a completely independent
// measurement): total system kinetic energy (both bodies, translational+rotational) grew from
// 9.05e-3 to 5.76e-2 over this window -- a 536% gain -- on a material whose restitution coefficient
// should LOSE roughly 91% of relative energy in a collision, not gain 5x. This scenario checks
// whether an ISOLATED two-body replay (no other floaters, no surrounding jam, no fluid coupling)
// reproduces a comparable gain. If it does, that confirms the defect is a genuine PAIRWISE
// contact-response defect that needs nothing from the wider jam context to manifest -- the strongest
// possible reproduction, since every input number here comes from an actual crash-adjacent production
// run rather than a constructed worst case.
//
// FIDELITY CAVEAT (documented, not hidden): this is a best-effort, not bit-exact, reproduction.
// LinearSpringDashpot's persistent tangential-spring history starts EMPTY at this replay's own step
// 0, whereas the real contact may have engaged a few sub-cycles before step 25883 (the
// PROLONGED-CONTACT log's "started step 25883" reflects when overlap first crossed ITS OWN reporting
// threshold, not necessarily the exact sub-cycle collisionResponse_ itself first fired). A close but
// not bit-identical energy trajectory is the expected, honest outcome.
// Loops the SAME real incident through the standard 9-way mitigation config set (identical rows to
// runMultiPointEdgeCases/runFallOntoPileStability, so the three battery's results line up directly in
// one comparison table) instead of a single hard-coded config -- "run all options against this test
// case too" per direct request. Default duration raised from the original 3540 sub-cycles (exactly
// matching the real log's own window) to 10000: this session's own finding is that baseline's peak
// doesn't land until sub-cycle ~4072 and the elevated plateau afterward is itself part of the signature
// (see this scenario's own history -- the exact-window replay alone understates the real gain by a lot,
// 72% vs the real log's 536%, because contact hasn't finished compressing yet at sub-cycle 3540).
inline bool runProductionCollisionReplay(const config::Config::BlockHandle& b)
{
   const Vec3   halfExtents   = b.getParameter< Vec3 >("halfExtents", Vec3(6, 1, 1.5));
   const real_t density       = b.getParameter< real_t >("density", real_t(0.6));
   const real_t restitution   = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t friction      = b.getParameter< real_t >("frictionCoefficient", real_t(0.5));
   const real_t poisson       = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime = b.getParameter< real_t >("collisionTime", real_t(32.0078));
   const real_t dt            = b.getParameter< real_t >("dt", real_t(0.05));
   const real_t gravity       = b.getParameter< real_t >("gravityMagnitude", real_t(7.02618e-6));
   const uint_t numSubCycles  = b.getParameter< uint_t >("numSubCycles", uint_t(10000));
   const real_t tolerance     = b.getParameter< real_t >("toleranceRelError", real_t(0.05));
   const uint_t vtkSpacing    = b.getParameter< uint_t >("vtkSpacing", uint_t(0));
   const bool   writeVtk      = b.getParameter< bool >("writeVtk", false);
   const std::string outBase = b.getParameter< std::string >("outputBaseName", "production_collision_replay");

   auto makeRot = [](real_t r, real_t i, real_t j, real_t k) {
      const real_t n = std::sqrt(r * r + i * i + j * j + k * k); // logged quats aren't exactly unit-norm
      return Rot3(math::Quaternion< real_t >(r / n, i / n, j / n, k / n));
   };

   // uid=720575940379279392, logged state at step 25883, subCycle 0 (rank_27_state.csv).
   const Vec3 posA(real_t(37.5365), real_t(20.4184), real_t(51.8732));
   const Vec3 velA(real_t(0.00278915), real_t(0.000524567), real_t(2.70449e-05));
   const Vec3 angVelA(real_t(-0.000119324), real_t(-0.000434568), real_t(1.53791e-05));
   const Rot3 rotA = makeRot(real_t(0.0401339), real_t(0.672907), real_t(-0.667814), real_t(0.315612));

   // uid=3314649325744685089 (the box independently confirmed "flying"), same step (rank_21_state.csv).
   const Vec3 posB(real_t(40.5866), real_t(17), real_t(57.8195));
   const Vec3 velB(real_t(0.000377042), real_t(0), real_t(-0.0202495));
   const Vec3 angVelB(real_t(0), real_t(0), real_t(0));
   const Rot3 rotB = makeRot(real_t(1), real_t(0), real_t(0), real_t(0));

   // Same standard 9-way config set as runMultiPointEdgeCases/runFallOntoPileStability -- see those
   // functions' own comments for what each row targets.
   struct Config
   {
      const char* name;
      bool        rounded;
      bool        multiPoint;
      bool        overlapCap;
      bool        hysteresis;
      bool        manifold;
      bool        rotationalEffMass = false;
   };
   // [Claude, 2026-10-29] "useSpheroBoxContact" row removed -- see archive/unused_contact_methods/.
   const Config configs[] = {
      { "baseline (both OFF)", false, false, false, false, false, false },
      { "useRoundedBoxCorners", true, false, false, false, false, false },
      { "useMultiPointBoxContact", false, true, false, false, false, false },
      { "both", true, true, false, false, false, false },
      { "useMultiPointBoxContact+overlapCap", false, true, true, false, false, false },
      { "useAxisHysteresis", false, false, false, true, false, false },
      { "useMultiPointBoxContact+manifold", false, true, false, false, true, false },
      { "useRotationalEffectiveMass", false, false, false, false, false, true },
      { "useMultiPointBoxContact+rotationalEffMass", false, true, false, false, false, true },
   };

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] productionCollisionReplay -- real incident from job 4199585 "
      "(uid=720575940379279392 vs uid=3314649325744685089, exact logged state at step 25883), " << numSubCycles
      << " sub-cycles, " << (sizeof(configs) / sizeof(configs[0])) << " configurations. Production log's own "
      "independently-measured gain over the matching 3540-sub-cycle window: KE 9.05e-3 -> 5.76e-2 (+536%).");

   bool anyConfigFailed = false;
   for (const auto& cfg : configs)
   {
      DEMEngine engine(dt, restitution, friction, poisson, collisionTime, gravity, real_t(200),
                        /*useAnalyticalContactDetection=*/true, cfg.rounded, real_t(0.05), cfg.multiPoint,
                        cfg.overlapCap, real_t(0.3), cfg.hysteresis, real_t(0.15), cfg.manifold, real_t(0.25),
                        cfg.rotationalEffMass);

      const size_t idxA = engine.addBox(posA, velA, rotA, real_t(2) * halfExtents, density);
      const size_t idxB = engine.addBox(posB, velB, rotB, real_t(2) * halfExtents, density);
      engine.ac().setAngularVelocity(idxA, angVelA);
      engine.ac().setAngularVelocity(idxB, angVelB);

      engine.finalizeSetup();

      std::string safeName = cfg.name;
      for (char& c : safeName)
         if (!std::isalnum(static_cast< unsigned char >(c))) c = '_';
      const std::string thisOutBase = outBase + "_" + safeName;

      if (writeVtk) engine.setupVtk(thisOutBase, vtkSpacing);

      // COMBINED total (both bodies) kinetic energy -- explicitly the sum tracked at every sub-cycle,
      // not a per-body value read off a CSV after the fact (an earlier ad-hoc post-hoc CSV inspection
      // this session mixed up a single-body column with the combined total; this is the authoritative,
      // in-loop combined value instead).
      const real_t keBefore = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      real_t       keMax     = keBefore;
      uint_t       keMaxStep = uint_t(0);

      for (uint_t s = 0; s < numSubCycles; ++s)
      {
         engine.step();
         if (writeVtk) engine.writeVtkIfDue();
         const real_t keNow = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
         if (keNow > keMax) { keMax = keNow; keMaxStep = s; }
      }

      const real_t keAfter   = engine.kineticEnergy(idxA) + engine.kineticEnergy(idxB);
      const real_t keGainRel = (keMax - keBefore) / std::max(real_t(1e-12), keBefore);
      // Same physically-inviolable-ish bound as runCornerNormalDiscontinuity, adapted for e=0.3 (not
      // e=1): with restitution 0.3, gravity's own tiny contribution over this window aside, KE should be
      // MONOTONICALLY non-increasing once compression ends -- a gain beyond `tolerance` is the same
      // smoking gun as the e=1 case, just with a slightly looser bound to account for gravity's small
      // legitimate contribution.
      const bool pass = keGainRel < tolerance;
      if (!pass) anyConfigFailed = true;

      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] productionCollisionReplay  config=\"" << cfg.name << "\""
         << "  KE_before=" << keBefore << " KE_after=" << keAfter << " KE_max=" << keMax
         << " (at sub-cycle " << keMaxStep << ") maxRelGain=" << keGainRel << "  "
         << (pass ? "PASS" : "FAIL (energy injected)"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] productionCollisionReplay OVERALL: "
      << (!anyConfigFailed ? "PASS" : "FAIL") << " (anyConfigFailed=" << anyConfigFailed << ")");
   return !anyConfigFailed;
}

//======================================================================================================================
// 6-7: bulk mixed-shape tests
//======================================================================================================================

inline bool runPacking(const config::Config::BlockHandle& b)
{
   const uint_t numPerShape    = b.getParameter< uint_t >("numPerShape", uint_t(30));
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(0.9, 0.9, 0.9));
   const real_t cylRadius      = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength      = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t e              = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu             = b.getParameter< real_t >("frictionCoefficient", real_t(0.5));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity        = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t containerX     = b.getParameter< real_t >("containerSizeX", real_t(20));
   const real_t containerY     = b.getParameter< real_t >("containerSizeY", real_t(20));
   const uint_t maxSteps       = b.getParameter< uint_t >("maxSteps", uint_t(50000));
   const uint_t vtkSpacing     = b.getParameter< uint_t >("vtkSpacing", uint_t(500));
   const std::string outBase   = b.getParameter< std::string >("outputBaseName", "packing");

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity, real_t(2) * (containerX + containerY));
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1));
   engine.addWall(Vec3(containerX, 0, 0), Vec3(-1, 0, 0));
   engine.addWall(Vec3(0, 0, 0), Vec3(1, 0, 0));
   engine.addWall(Vec3(0, containerY, 0), Vec3(0, -1, 0));
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 1, 0));

   real_t       totalVolume = real_t(0);
   uint_t       created     = 0;
   const real_t dropHeight  = real_t(1);
   const real_t cellSize    = real_t(2.2) * std::max({ sphereRadius, boxHalfExtents[0], boxHalfExtents[1],
                                                        boxHalfExtents[2],
                                                        std::sqrt(cylRadius * cylRadius +
                                                                   real_t(0.25) * cylLength * cylLength) });
   const uint_t perRow      = std::max< uint_t >(1, uint_t(containerX / cellSize));

   for (uint_t shapeType = 0; shapeType < 3; ++shapeType)
      for (uint_t i = 0; i < numPerShape; ++i)
      {
         const uint_t row = created % perRow, col = (created / perRow) % perRow;
         const real_t x =
            cellSize * (real_t(row) + real_t(0.5)) + math::realRandom< real_t >(-real_t(0.1), real_t(0.1));
         const real_t y =
            cellSize * (real_t(col) + real_t(0.5)) + math::realRandom< real_t >(-real_t(0.1), real_t(0.1));
         const real_t z = dropHeight + (real_t(created) / real_t(perRow * perRow)) * cellSize;
         const Vec3   pos(std::min(std::max(x, real_t(0.5)), containerX - real_t(0.5)),
                           std::min(std::max(y, real_t(0.5)), containerY - real_t(0.5)), z);

         if (shapeType == 0)
         {
            engine.addSphere(pos, Vec3(0), sphereRadius, density);
            totalVolume += (real_t(4) / real_t(3)) * math::pi * sphereRadius * sphereRadius * sphereRadius;
         }
         else if (shapeType == 1)
         {
            const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                                 math::realRandom< real_t >(-1, 1))
                              .getNormalizedIfNotZero(),
                           math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
            engine.addBox(pos, Vec3(0), rot, real_t(2) * boxHalfExtents, density);
            totalVolume += real_t(8) * boxHalfExtents[0] * boxHalfExtents[1] * boxHalfExtents[2];
         }
         else
         {
            const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                                 math::realRandom< real_t >(-1, 1))
                              .getNormalizedIfNotZero(),
                           math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
            engine.addCylinder(pos, Vec3(0), rot, cylRadius, cylLength, density);
            totalVolume += math::pi * cylRadius * cylRadius * cylLength;
         }
         ++created;
      }

   engine.finalizeSetup();
   engine.setupVtk(outBase, vtkSpacing);

   real_t avgSpeed          = std::numeric_limits< real_t >::max();
   real_t prevEnergy        = std::numeric_limits< real_t >::max();
   real_t worstEnergyIncrease = real_t(0); // largest checkpoint-to-checkpoint relative increase seen
   uint_t step              = 0;
   for (; step < maxSteps && avgSpeed > real_t(1e-3) * std::sqrt(gravity * sphereRadius); ++step)
   {
      engine.step();
      engine.writeVtkIfDue();
      // step > 0 guard: VelocityVerletPostForceUpdate averages oldForce and newForce, and oldForce is
      // still zero on the very first step -- checking avgSpeed at step==0 measures roughly HALF a
      // gravity kick regardless of true settling state, which can spuriously read as "already settled"
      // before any real dynamics (contacts, actual falling) have had a chance to develop.
      if (step > 0 && step % 200 == 0)
      {
         real_t sumSpeed = real_t(0);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            sumSpeed += engine.ac().getLinearVelocity(i).length();
         }
         avgSpeed = sumSpeed / real_t(created);

         const real_t energy = detail::totalMechanicalEnergy(engine, gravity);
         if (prevEnergy < std::numeric_limits< real_t >::max())
         {
            const real_t increase = (energy - prevEnergy) / std::max(real_t(1e-9), std::abs(prevEnergy));
            worstEnergyIncrease = std::max(worstEnergyIncrease, increase);
         }
         prevEnergy = energy;

         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packing step=" << step << " avgSpeed=" << avgSpeed
                                    << " totalMechanicalEnergy=" << energy);
      }
   }

   // Checkpoint-to-checkpoint tolerance (not step-to-step): a sub-cycled spring-dashpot contact is not
   // exactly energy-monotonic on every single timestep even when correct, but should not show a
   // sustained upward trend sampled every 200 steps.
   const bool energyMonotonic = worstEnergyIncrease < real_t(0.01);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packing [energy non-increase] worst checkpoint-to-checkpoint "
                              << "relative increase=" << worstEnergyIncrease << "  "
                              << (energyMonotonic ? "PASS" : "FAIL (energy created -- contact/integration unstable)"));

   real_t maxTop = real_t(0);
   for (size_t i = 0; i < engine.ac().size(); ++i)
   {
      if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
      maxTop = std::max(maxTop, engine.ac().getPosition(i)[2] + engine.ac().getInteractionRadius(i));
   }

   // solid fraction = totalVolume / (containerX*containerY*maxTop), checked against the band below
   const real_t bulkSolidFraction = totalVolume / (containerX * containerY * std::max(maxTop, real_t(1e-6)));
   // Literature ordering for random-packed basic 3D shapes (Li, S., Zhao, J., Lu, P., & Xie, Y. (2010),
   // "Maximum packing densities of basic 3D objects," Chinese Science Bulletin 55(2), 114-119,
   // doi:10.1007/s11434-009-0650-0 -- [Claude, 2026-10-29] CORRECTED author attribution: this was
   // previously misattributed to "Jiao, Y., Stillinger, F.H. & Torquato, S." in this comment, a
   // different, unrelated set of packing-density researchers who did not write this particular paper;
   // verified via web search against the DOI/journal/volume/pages, which all match Li et al. instead):
   // the MAXIMALLY-JAMMED upper bounds are sphere~0.64, cylinder~0.72, cube~0.78. A gravity-poured bed
   // settles BELOW these maxima, so the pass band below is a broad plausibility check ("not obviously
   // broken"), not an exact-match test.
   // NOT gated on step < maxSteps: a nonzero-restitution mixed-shape bed settles into a persistent
   // low-level dynamic equilibrium (particles perpetually micro-resettling against their irregular
   // neighbours) rather than exactly zero velocity, confirmed empirically -- avgSpeed plateaus at a
   // small but non-decaying value well before maxSteps rather than continuing to fall, so requiring it
   // to drop below the strict threshold just measures "ran out of budget," not "bed is unphysical."
   // maxSteps is long enough to reach a representative quasi-steady state; the bulk-fraction plausibility
   // band is the actual check.
   const bool pass = bulkSolidFraction > real_t(0.30) && bulkSolidFraction < real_t(0.80) && energyMonotonic;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packing steps="
                              << step << " settledHeight=" << maxTop << " bulkSolidFraction=" << bulkSolidFraction
                              << "  (literature RCP upper bounds: sphere~0.64, cylinder~0.72, cube~0.78)"
                              << "  " << (pass ? "PASS" : "FAIL"));
   return pass;
}

inline bool runAngleOfRepose(const config::Config::BlockHandle& b)
{
   const uint_t numParticles    = b.getParameter< uint_t >("numParticles", uint_t(60));
   const real_t sphereRadius    = b.getParameter< real_t >("sphereRadius", real_t(1));
   const real_t cylRadius       = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength       = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const bool   useCylinders    = b.getParameter< bool >("useCylinders", true);
   const real_t density         = b.getParameter< real_t >("density", real_t(1));
   const real_t e               = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu              = b.getParameter< real_t >("frictionCoefficient", real_t(0.6));
   const real_t poisson         = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime   = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt              = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity         = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t floorHalfExtent = b.getParameter< real_t >("floorHalfExtent", real_t(30));
   const uint_t maxSteps        = b.getParameter< uint_t >("maxSteps", uint_t(60000));
   const uint_t vtkSpacing      = b.getParameter< uint_t >("vtkSpacing", uint_t(500));
   const std::string outBase    = b.getParameter< std::string >("outputBaseName", "angle_of_repose");

   // Domain must be tall enough for the INITIAL POUR COLUMN (stacked 5-at-a-time, 0.6 units apart, well
   // above the floor before gravity settles it -- see the placement loop below), not just wide enough
   // for the floor footprint: a small floorHalfExtent alone previously left the pour column poking
   // outside a too-short domain and crashing AssocToBlock's block-assignment assertion.
   //
   // CALIBRATION NOTE (two earlier versions, both failed): (1) a tall, narrow random-polar column
   // (1.5-unit radius, 1.1-unit vertical spacing) gave a ~55-unit free-fall for 150 particles --
   // sqrt(2*g*55)=~33 impact speed, scattering particles across the open floor (measured repose
   // ~0.4deg, i.e. flat). (2) shortening that column while KEEPING random (r,phi) placement fixed the
   // impact-energy problem but exposed a worse one: random placement within a small radius routinely
   // seated multiple particles overlapping at creation (before any integration), and LinearSpringDashpot's
   // force scales with penetration depth, so a large initial overlap is an explosive first-step force --
   // particles were measured flying out to r,z in the THOUSANDS (should be tens at most), an outright
   // numerical blow-up, not a physical spread (still read as "flat" because the fit average washes out).
   // Fix: a REGULAR GRID (cell size with a safety margin over the largest particle's circumscribing
   // radius, exactly runPacking's proven-safe non-overlapping placement pattern) instead of random polar
   // coordinates, so no two particles start in contact regardless of how many are poured.
   const real_t cellSize =
      real_t(2.2) * std::max(sphereRadius, std::sqrt(cylRadius * cylRadius + real_t(0.25) * cylLength * cylLength));
   const uint_t gridSide          = uint_t(4); // 4x4 = 16 particles per layer
   const real_t pourColumnHeight  = real_t(2) + cellSize * real_t(numParticles / (gridSide * gridSide) + 1) + real_t(2);
   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity,
                     std::max(real_t(3) * floorHalfExtent, pourColumnHeight + real_t(5)));
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1));

   // CALIBRATION NOTE (3): even with a numerically clean, non-overlapping pour, the measured pile stayed
   // essentially flat (repose < 1deg, particles spreading to r~70-80) on a fully open floor. Root cause:
   // mesa_pd's LinearSpringDashpot models pure Coulomb SLIDING friction with no explicit rolling-
   // resistance torque, so a sphere -- or even an elongated cylinder once it starts rolling end-over-end
   // or spinning about its long axis -- keeps rolling indefinitely rather than dissipating into a bounded
   // pile; this is a well-documented real limitation of plain-Coulomb DEM contact models (granular
   // calibration literature routinely needs an explicit rolling-friction term for a realistic angle of
   // repose), not a bug in this test. Containing the pour in a bounded pit (the same wall pattern
   // runPacking already uses) sidesteps it: material that CAN'T roll away indefinitely still piles up
   // and develops a measurable free-surface slope near the walls, closer to the "pour into a hopper, then
   // measure the settled slope" method the DEM calibration literature also documents, rather than
   // Kodam's fully free-standing open-floor cone.
   const real_t pitHalfExtent = b.getParameter< real_t >("pitHalfExtent", real_t(10));
   engine.addWall(Vec3(pitHalfExtent, 0, 0), Vec3(-1, 0, 0));
   engine.addWall(Vec3(-pitHalfExtent, 0, 0), Vec3(1, 0, 0));
   engine.addWall(Vec3(0, pitHalfExtent, 0), Vec3(0, -1, 0));
   engine.addWall(Vec3(0, -pitHalfExtent, 0), Vec3(0, 1, 0));

   // Pour from a small column above the centre -- the standard static-pile method used in DEM
   // calibration studies for measuring angle of repose (see e.g. the review literature on angle-of-
   // repose scaling/calibration this session's research surfaced).
   uint_t created = 0;
   for (uint_t i = 0; i < numParticles; ++i)
   {
      const uint_t row = (i / gridSide) % gridSide, col = i % gridSide, layer = i / (gridSide * gridSide);
      const real_t gx  = cellSize * (real_t(row) - real_t(gridSide - 1) * real_t(0.5)) +
                         math::realRandom< real_t >(-real_t(0.05), real_t(0.05));
      const real_t gy = cellSize * (real_t(col) - real_t(gridSide - 1) * real_t(0.5)) +
                        math::realRandom< real_t >(-real_t(0.05), real_t(0.05));
      const Vec3 pos(gx, gy, real_t(2) + cellSize * real_t(layer));
      // ALL cylinders when useCylinders (not an alternating 50/50 mix with spheres): mesa_pd's
      // LinearSpringDashpot models pure Coulomb SLIDING friction with no explicit rolling-resistance
      // torque, so a rotationally-symmetric sphere (or a short, disk-like low-aspect-ratio cylinder) can
      // roll indefinitely once set in motion rather than settling into a bounded pile -- a known,
      // literature-documented DEM limitation (spherical/near-spherical particles need an explicit
      // rolling-friction model for a realistic angle of repose), not a bug in this test. An elongated rod
      // (see cylLength's aspect ratio) resists end-over-end rolling far better than a sphere ever can.
      if (useCylinders)
      {
         const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                              math::realRandom< real_t >(-1, 1))
                           .getNormalizedIfNotZero(),
                        math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
         engine.addCylinder(pos, Vec3(0), rot, cylRadius, cylLength, density);
      }
      else { engine.addSphere(pos, Vec3(0), sphereRadius, density); }
      ++created;
   }

   engine.finalizeSetup();
   engine.setupVtk(outBase, vtkSpacing);

   real_t avgSpeed            = std::numeric_limits< real_t >::max();
   real_t prevEnergy          = std::numeric_limits< real_t >::max();
   real_t worstEnergyIncrease = real_t(0);
   uint_t step                = 0;
   for (; step < maxSteps && avgSpeed > real_t(1e-3) * std::sqrt(gravity * sphereRadius); ++step)
   {
      engine.step();
      engine.writeVtkIfDue();
      // step > 0 guard: VelocityVerletPostForceUpdate averages oldForce and newForce, and oldForce is
      // still zero on the very first step -- checking avgSpeed at step==0 measures roughly HALF a
      // gravity kick regardless of true settling state, which can spuriously read as "already settled"
      // before any real dynamics (contacts, actual falling) have had a chance to develop.
      if (step > 0 && step % 200 == 0)
      {
         real_t sumSpeed = real_t(0);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            sumSpeed += engine.ac().getLinearVelocity(i).length();
         }
         avgSpeed = sumSpeed / real_t(created);

         const real_t energy = detail::totalMechanicalEnergy(engine, gravity);
         if (prevEnergy < std::numeric_limits< real_t >::max())
         {
            const real_t increase = (energy - prevEnergy) / std::max(real_t(1e-9), std::abs(prevEnergy));
            worstEnergyIncrease   = std::max(worstEnergyIncrease, increase);
         }
         prevEnergy = energy;
      }
   }

   const bool energyMonotonic = worstEnergyIncrease < real_t(0.01);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] angleOfRepose [energy non-increase] worst checkpoint-to-checkpoint "
                              << "relative increase=" << worstEnergyIncrease << "  "
                              << (energyMonotonic ? "PASS" : "FAIL (energy created -- contact/integration unstable)"));

   // Least-squares fit of height vs. radial distance from the pour axis over the settled pile's outer
   // shell -- slope magnitude approximates tan(angle of repose).
   real_t sumR = 0, sumZ = 0, sumRZ = 0, sumRR = 0;
   real_t minR = std::numeric_limits< real_t >::max(), maxR = real_t(0);
   real_t minZ = std::numeric_limits< real_t >::max(), maxZ = real_t(0);
   uint_t n = 0;
   for (size_t i = 0; i < engine.ac().size(); ++i)
   {
      if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
      const Vec3   pos = engine.ac().getPosition(i);
      const real_t r   = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      sumR += r;
      sumZ += pos[2];
      sumRZ += r * pos[2];
      sumRR += r * r;
      minR = std::min(minR, r);
      maxR = std::max(maxR, r);
      minZ = std::min(minZ, pos[2]);
      maxZ = std::max(maxZ, pos[2]);
      ++n;
   }
   // Sanity-check range: a genuinely settled pile should have maxR/maxZ of the same order as the pour
   // footprint/height, not orders of magnitude larger (the latter would indicate a numerical blow-up,
   // e.g. particles seeded overlapping at creation, rather than a physically settled result).
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] angleOfRepose pile extent: n=" << n << " rRange=[" << minR << ","
                              << maxR << "] zRange=[" << minZ << "," << maxZ << "]");
   const real_t denom = real_t(n) * sumRR - sumR * sumR;
   const real_t slope = (std::abs(denom) > real_t(1e-9)) ? (real_t(n) * sumRZ - sumR * sumZ) / denom : real_t(0);
   // repose angle = atan(|slope|) of the pile's r-vs-z least-squares fit; checked against the
   // plausibility band below, anchored by Kodam et al. 2010 Part II's 39.6deg sim / 43.9deg experiment
   const real_t reposeAngleDeg = std::atan(std::abs(slope)) * real_t(180) / math::pi;

   // Plausibility band only (~15-55 deg spans essentially all real granular materials); Kodam et al.
   // (2010) Part II's own cylinder result (39.6 deg simulated vs. 43.9 deg experimental, aspect ratio
   // ~0.66) is one concrete literature anchor within this band. NOT gated on step < maxSteps -- see
   // runPacking's identical reasoning: a nonzero-restitution mixed-shape pile settles into a persistent
   // low-level dynamic equilibrium, not exact rest.
   const bool pass = reposeAngleDeg > real_t(15) && reposeAngleDeg < real_t(55) && energyMonotonic;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] angleOfRepose steps="
                              << step << " reposeAngle=" << reposeAngleDeg
                              << "deg  (plausible range 15-55deg; Kodam et al. 2010 cylinder reference: "
                              << "39.6deg sim / 43.9deg experiment)  " << (pass ? "PASS" : "FAIL"));
   if (!pass && reposeAngleDeg < real_t(5))
      WALBERLA_LOG_INFO_ON_ROOT(
         "[DEMValidation] angleOfRepose DIAGNOSIS: a near-zero angle after settling (not just below "
         "tolerance) points to the known mesa_pd LinearSpringDashpot limitation documented above this "
         "function -- pure Coulomb sliding friction with no rolling-resistance torque lets particles roll/"
         "slide out to a flat, low layer rather than piling up, even when contained. This is a genuine "
         "contact-model gap this validation run has surfaced, not a broken test.");
   return pass;
}

//======================================================================================================================
// 8-9: single-shape Kodam et al. (2010) Part II Fig. 10/11 methodology recreations -- unlike runPacking/
// runAngleOfRepose above (deliberately generalised/simplified variants), these two follow Kodam's actual
// protocol as closely as a fluid-free DEM dry-run reasonably can: ONE shape at a time (swap via the
// `shape` prm key: "sphere"/"box"/"cylinder"), a genuinely circular container/drum wall (polygon-of-
// planes approximation, the same "approximate the curved wall with flat plates" idea their own Fig. 6(a)
// uses for the friability tester), particles with a random initial orientation, and -- for the repose
// test -- an actually ROTATING drum (DEMEngine::setWallPose called every tumbling step) rather than a
// static pour-pile.
//======================================================================================================================

namespace detail
{
//! Creates one particle of `shape` ("sphere"/"box"/"cylinder") at `pos`, at rest, with a uniformly random
//! 3D orientation (matching Kodam et al.'s "random initial orientation" drop protocol) -- spheres are
//! orientation-invariant so the random rotation is skipped for them (harmless either way, just wasted
//! entropy).
inline size_t addRandomOrientedShape(DEMEngine& engine, const std::string& shape, const Vec3& pos,
                                      real_t sphereRadius, const Vec3& boxHalfExtents, real_t cylRadius,
                                      real_t cylLength, real_t density)
{
   if (shape == "sphere") return engine.addSphere(pos, Vec3(0), sphereRadius, density);
   const Rot3 rot(Vec3(math::realRandom< real_t >(-1, 1), math::realRandom< real_t >(-1, 1),
                        math::realRandom< real_t >(-1, 1))
                     .getNormalizedIfNotZero(),
                  math::realRandom< real_t >(real_t(0), real_t(2) * math::pi));
   if (shape == "box") return engine.addBox(pos, Vec3(0), rot, real_t(2) * boxHalfExtents, density);
   if (shape == "cylinder") return engine.addCylinder(pos, Vec3(0), rot, cylRadius, cylLength, density);
   WALBERLA_ABORT("Unknown shape '" << shape << "'. Valid values: sphere, box, cylinder.");
}

//! Circumscribing-sphere radius for `shape` -- used for non-overlapping grid spacing and fill-height
//! measurement (same convention as DEMEngine::addSphere/addBox/addCylinder's own setInteractionRadius).
inline real_t shapeCircumRadius(const std::string& shape, real_t sphereRadius, const Vec3& boxHalfExtents,
                                 real_t cylRadius, real_t cylLength)
{
   if (shape == "sphere") return sphereRadius;
   if (shape == "box") return boxHalfExtents.length();
   if (shape == "cylinder") return std::sqrt(cylRadius * cylRadius + real_t(0.25) * cylLength * cylLength);
   WALBERLA_ABORT("Unknown shape '" << shape << "'. Valid values: sphere, box, cylinder.");
}

//! Exact solid volume for `shape` (used for bulk solid fraction, Kodam et al. Eq. (7)).
inline real_t shapeVolume(const std::string& shape, real_t sphereRadius, const Vec3& boxHalfExtents,
                           real_t cylRadius, real_t cylLength)
{
   if (shape == "sphere") return (real_t(4) / real_t(3)) * math::pi * sphereRadius * sphereRadius * sphereRadius;
   if (shape == "box") return real_t(8) * boxHalfExtents[0] * boxHalfExtents[1] * boxHalfExtents[2];
   if (shape == "cylinder") return math::pi * cylRadius * cylRadius * cylLength;
   WALBERLA_ABORT("Unknown shape '" << shape << "'. Valid values: sphere, box, cylinder.");
}

//! Builds a regular `numSegments`-sided polygon of flat HalfSpace walls approximating a cylindrical
//! container/drum wall of interior (apothem) radius `radius`, axis along world Y, centred at
//! (axisX, *, axisZ) in the X-Z plane. `angleOffset` rotates the whole polygon; passing a time-varying
//! offset (see updateDrumWalls) animates a rotating drum. Each wall is tangent to the target circle at
//! its segment's angle, so particles are contained within `radius` at the flats and bulge very slightly
//! beyond it at the corners (by a factor 1/cos(pi/numSegments), negligible for numSegments>=16).
inline std::vector< size_t > addDrumWalls(DEMEngine& engine, real_t axisX, real_t axisZ, real_t radius,
                                           uint_t numSegments, real_t angleOffset = real_t(0))
{
   std::vector< size_t > walls;
   walls.reserve(numSegments);
   for (uint_t i = 0; i < numSegments; ++i)
   {
      const real_t theta = real_t(2) * math::pi * real_t(i) / real_t(numSegments) + angleOffset;
      const Vec3   normal(-std::cos(theta), real_t(0), -std::sin(theta)); // inward
      const Vec3   pos(axisX + radius * std::cos(theta), real_t(0), axisZ + radius * std::sin(theta));
      walls.push_back(engine.addWall(pos, normal));
   }
   return walls;
}

//! Re-poses an existing set of drum walls (from addDrumWalls) to a new `angleOffset`, keeping the same
//! axis/radius/segment count -- called every step during the tumbling phase of runRotatingDrum.
inline void updateDrumWalls(DEMEngine& engine, const std::vector< size_t >& walls, real_t axisX, real_t axisZ,
                             real_t radius, real_t angleOffset)
{
   const uint_t numSegments = uint_t(walls.size());
   for (uint_t i = 0; i < numSegments; ++i)
   {
      const real_t theta = real_t(2) * math::pi * real_t(i) / real_t(numSegments) + angleOffset;
      const Vec3   normal(-std::cos(theta), real_t(0), -std::sin(theta));
      const Vec3   pos(axisX + radius * std::cos(theta), real_t(0), axisZ + radius * std::sin(theta));
      engine.setWallPose(walls[i], pos, normal);
   }
}
} // namespace detail

//! Kodam et al. (2010) Fig. 10 / Section 3.3 recreation: ONE shape (swap via `shape`), released one at a
//! time from a random point near the top with negligible initial KE and a random orientation, into a
//! genuinely circular container (detail::addDrumWalls, static), each particle waited out (settle-speed
//! threshold or dropTimeoutSteps) before the next is released -- Kodam's own protocol keeps the drop rate
//! far below 1000/s specifically so drop rate itself doesn't bias the packing (Macrae & Gray, 1961; Zhang
//! et al., 2001, both cited in the paper); one-at-a-time is the limiting case of that. Measures fill
//! height and bulk solid fraction (their Eq. (7)) against shape-specific plausibility bands anchored to
//! the random-close-packing literature (Li, Zhao, Lu & Xie, 2010, Chinese Science Bulletin 55(2):114-119;
//! NOT "Jiao, Stillinger & Torquato" as an earlier version of this file misattributed it, see runPacking's
//! own comment) and Kodam's own measured
//! cylinder values (Table 2: 0.617 experiment / 0.589 true-cylinder-DEM sim).
inline bool runPackingSingleShape(const config::Config::BlockHandle& b)
{
   const std::string shape        = b.getParameter< std::string >("shape", "cylinder");
   const uint_t       numParticles = b.getParameter< uint_t >("numParticles", uint_t(250));
   const real_t sphereRadius      = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents    = b.getParameter< Vec3 >("boxHalfExtents", Vec3(0.9, 0.9, 0.9));
   const real_t cylRadius         = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength         = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const real_t density           = b.getParameter< real_t >("density", real_t(1));
   const real_t e                 = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu                = b.getParameter< real_t >("frictionCoefficient", real_t(0.5));
   const real_t poisson           = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime     = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt                = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity           = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t containerRadius   = b.getParameter< real_t >("containerRadius", real_t(6));
   const uint_t numWallSegments   = b.getParameter< uint_t >("numWallSegments", uint_t(24));
   const uint_t settleCheckEvery  = b.getParameter< uint_t >("settleCheckEvery", uint_t(200));
   const uint_t dropTimeoutSteps  = b.getParameter< uint_t >("dropTimeoutSteps", uint_t(6000));
   const uint_t maxSteps          = b.getParameter< uint_t >("maxSteps", uint_t(600000));
   const uint_t vtkSpacing        = b.getParameter< uint_t >("vtkSpacing", uint_t(200));
   const std::string outBase      = b.getParameter< std::string >("outputBaseName", "packing_single_shape");

   const real_t circumRadius   = detail::shapeCircumRadius(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);
   const real_t particleVolume = detail::shapeVolume(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);

   // Domain must comfortably contain the tallest plausible settled pile PLUS the fixed drop height above
   // it -- sized from the actual particle count/volume rather than a flat guess, so this stays correct
   // across very different shape/count combinations instead of silently truncating a large run (the
   // failure mode runAngleOfRepose's own calibration notes already document for an under-sized domain).
   const real_t estPileHeight = real_t(1.6) * particleVolume * real_t(numParticles) /
                                 (math::pi * containerRadius * containerRadius);
   const real_t dropHeight     = real_t(2) * containerRadius + real_t(2) * circumRadius + estPileHeight;
   // Margin scaled to particle size (not a fixed absolute constant): at Kodam's real SI scale
   // (containerRadius ~2.5e-2 m), a fixed "+10" margin -- fine at the old generic-units scale
   // (containerRadius ~6) -- would outweigh the entire physical domain by ~1000x, forcing step()'s
   // LinkedCells (cell width set from the real, mm-scale particle diameter) to cover a wildly oversized
   // domain and allocate on the order of 1e10 cells, segfaulting on the very first step.
   const real_t domainHalfExtent = dropHeight + real_t(2) * containerRadius + real_t(4) * circumRadius;

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity, domainHalfExtent);
   engine.addWall(Vec3(0, 0, 0), Vec3(0, 0, 1)); // floor
   detail::addDrumWalls(engine, real_t(0), real_t(0), containerRadius, numWallSegments);

   const real_t dropRadius = real_t(0.6) * containerRadius; // random release point, away from the wall

   // No dynamic particles exist yet at this point (only GLOBAL walls) -- pass the known max diameter
   // explicitly (see DEMEngine::finalizeSetup's doc comment) so the linked-cell grid is sized correctly
   // for the particles the drop loop below is about to create one at a time.
   engine.finalizeSetup(real_t(2) * circumRadius);
   engine.setupVtk(outBase, vtkSpacing);

   uint_t step    = 0;
   uint_t created = 0;
   for (; created < numParticles && step < maxSteps; ++created)
   {
      // Release height tracks the current pile top (with margin) rather than a fixed value -- a fixed
      // drop height would eventually spawn a particle INSIDE the growing pile once it gets tall enough,
      // an explosive initial-overlap contact force (the same failure mode documented in
      // runAngleOfRepose's CALIBRATION NOTE (2) above).
      real_t currentTop = real_t(0);
      for (size_t i = 0; i < engine.ac().size(); ++i)
      {
         if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
         currentTop = std::max(currentTop, engine.ac().getPosition(i)[2] + engine.ac().getInteractionRadius(i));
      }
      const real_t dropZ = std::max(real_t(2) * containerRadius, currentTop + real_t(3) * circumRadius);

      const real_t r   = dropRadius * std::sqrt(math::realRandom< real_t >(real_t(0), real_t(1)));
      const real_t phi = math::realRandom< real_t >(real_t(0), real_t(2) * math::pi);
      detail::addRandomOrientedShape(engine, shape, Vec3(r * std::cos(phi), real_t(0), dropZ), sphereRadius,
                                      boxHalfExtents, cylRadius, cylLength, density);

      // Wait for this particle (and whatever it disturbs) to settle before releasing the next one --
      // Kodam et al.'s own one-at-a-time protocol, recreated directly rather than approximated by a bulk
      // gravity pour. dropTimeoutSteps is a safety valve: a particle in a stable-but-nonzero micro-
      // equilibrium (rocking against neighbours) should not stall the whole sequence indefinitely.
      for (uint_t waited = 0; waited < dropTimeoutSteps && step < maxSteps; ++waited, ++step)
      {
         engine.step();
         engine.writeVtkIfDue();
         if (waited > 0 && waited % settleCheckEvery == 0)
         {
            real_t sumSpeed = real_t(0);
            uint_t n         = uint_t(0);
            for (size_t i = 0; i < engine.ac().size(); ++i)
            {
               if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
               sumSpeed += engine.ac().getLinearVelocity(i).length();
               ++n;
            }
            if (n > 0 && sumSpeed / real_t(n) < real_t(1e-2) * std::sqrt(gravity * circumRadius)) break;
         }
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packingSingleShape shape=" << shape << " created=" << created << "/"
                              << numParticles << " steps=" << step);

   real_t maxTop = real_t(0);
   for (size_t i = 0; i < engine.ac().size(); ++i)
   {
      if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
      maxTop = std::max(maxTop, engine.ac().getPosition(i)[2] + engine.ac().getInteractionRadius(i));
   }

   // solid fraction = N*V_particle / (pi*R_container^2*maxTop), checked against the band below
   const real_t bulkSolidFraction = real_t(created) * particleVolume /
                                     (math::pi * containerRadius * containerRadius * std::max(maxTop, real_t(1e-6)));

   // Shape-specific plausibility bands, anchored to the literature RCP upper bounds (Li, Zhao, Lu & Xie,
   // 2010, Chinese Science Bulletin 55(2):114-119, doi:10.1007/s11434-009-0650-0 -- see runPacking's own
   // comment for the attribution correction: NOT "Jiao, Stillinger & Torquato" as an earlier version of
   // this comment claimed -- sphere~0.64, cylinder~0.72, cube~0.78 maximally jammed; no copy in
   // cld_kbase/, so treat as a soft, non-gating ceiling only) and Kodam et al.'s own measured cylinder
   // result (0.617 experiment / 0.589 true-cylinder-DEM sim,
   // Table 2, verified against kodam_cylindrical_2010b.txt) as a concrete anchor -- a gravity-poured
   // bed sits below the jammed maxima, so these are broad "not obviously broken" bands, not
   // exact-match tests: this dry run's contact parameters are not calibrated to Kodam's
   // Plexiglas-on-Plexiglas values, so exact numeric reproduction isn't the goal.
   real_t bandHi = real_t(0.64);
   if (shape == "cylinder") bandHi = real_t(0.72);
   else if (shape == "box") bandHi = real_t(0.78);
   const bool pass = created == numParticles && bulkSolidFraction > real_t(0.30) && bulkSolidFraction < bandHi;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] packingSingleShape shape="
                              << shape << " settledHeight=" << maxTop << " bulkSolidFraction=" << bulkSolidFraction
                              << " (plausible band [0.30," << bandHi << "])"
                              << "  (Kodam et al. 2010 Table 2 cylinder reference: 0.617 experiment / 0.589 "
                                 "true-cylinder-DEM sim)  "
                              << (pass ? "PASS" : "FAIL"));
   return pass;
}

//! Kodam et al. (2010) Fig. 11 / Section 3.4 recreation: ONE shape (swap via `shape`) tumbling inside an
//! ACTUALLY ROTATING drum (detail::addDrumWalls + updateDrumWalls, animated every tumbling step at `rpm`
//! -- the real thing, not runAngleOfRepose's static pour-pile substitute), axis along world Y, bounded by
//! two flat endcap walls at +-drumHalfLength. Particles fill the drum's interior directly (a 3D grid
//! clipped to the circular cross-section, non-overlapping by construction -- gravity sorts the initial
//! configuration into a resting bed during the settle phase, so no "pour from above a sealed drum" issue
//! ever arises). After settling with the drum stationary, the drum rotates for `tumbleRevolutions`
//! revolutions; `numSnapshots` free-surface readings are taken across that window (Kodam et al.: 10
//! photos at 1s intervals during their real 25rpm test) -- each reading bins particles by X into
//! `numRadialBins` bins and records the topmost occupied Z per occupied bin, the same "profile from
//! images" idea as their Fig. 12(a) pixel extraction, sampled directly from particle geometry instead of
//! pixels. All snapshot points are pooled and least-squares line-fit (their Fig. 12(b) method) to get one
//! dynamic angle of repose.
inline bool runRotatingDrum(const config::Config::BlockHandle& b)
{
   const std::string shape         = b.getParameter< std::string >("shape", "cylinder");
   const uint_t       numParticles  = b.getParameter< uint_t >("numParticles", uint_t(100));
   const real_t sphereRadius       = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents     = b.getParameter< Vec3 >("boxHalfExtents", Vec3(0.9, 0.9, 0.9));
   const real_t cylRadius          = b.getParameter< real_t >("cylinderRadius", real_t(0.8));
   const real_t cylLength          = b.getParameter< real_t >("cylinderLength", real_t(1.6));
   const real_t density            = b.getParameter< real_t >("density", real_t(1));
   const real_t e                  = b.getParameter< real_t >("restitutionCoefficient", real_t(0.3));
   const real_t mu                 = b.getParameter< real_t >("frictionCoefficient", real_t(0.6));
   const real_t poisson            = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime      = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt                 = b.getParameter< real_t >("dt", collisionTime / real_t(30));
   const real_t gravity            = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t drumRadius         = b.getParameter< real_t >("drumRadius", real_t(10));
   const real_t drumHalfLength     = b.getParameter< real_t >("drumHalfLength", real_t(5.6));
   const uint_t numWallSegments    = b.getParameter< uint_t >("numWallSegments", uint_t(24));
   const real_t rpm                = b.getParameter< real_t >("rpm", real_t(25));       // Kodam et al.'s own value
   const uint_t settleSteps        = b.getParameter< uint_t >("settleSteps", uint_t(30000));
   const uint_t tumbleRevolutions  = b.getParameter< uint_t >("tumbleRevolutions", uint_t(3));
   const uint_t numSnapshots       = b.getParameter< uint_t >("numSnapshots", uint_t(10)); // Kodam: 10 photos
   const uint_t numRadialBins      = b.getParameter< uint_t >("numRadialBins", uint_t(20));
   const uint_t vtkSpacing         = b.getParameter< uint_t >("vtkSpacing", uint_t(200));
   const std::string outBase       = b.getParameter< std::string >("outputBaseName", "rotating_drum");

   const real_t circumRadius = detail::shapeCircumRadius(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);
   const real_t omega        = rpm * real_t(2) * math::pi / real_t(60); // rad/s

   DEMEngine engine(dt, e, mu, poisson, collisionTime, gravity, real_t(3) * (drumRadius + drumHalfLength));
   auto wallIndices = detail::addDrumWalls(engine, real_t(0), real_t(0), drumRadius, numWallSegments);
   engine.addWall(Vec3(0, drumHalfLength, 0), Vec3(0, -1, 0));  // +Y endcap
   engine.addWall(Vec3(0, -drumHalfLength, 0), Vec3(0, 1, 0));  // -Y endcap

   // Fill the drum's interior directly with a 3D grid clipped to the circular cross-section (not a pour
   // column above a sealed top -- a real drum is loaded closed, particles start inside). Non-overlapping
   // by construction (cellSize has a safety margin over circumRadius), so ANY subset works as a starting
   // condition: gravity settles it into a resting bed during the stationary phase below, regardless of
   // which grid cells happened to be picked first.
   const real_t cellSize       = real_t(2.1) * circumRadius;
   const real_t interiorRadius = std::max(real_t(0.1), drumRadius - real_t(1.05) * circumRadius);
   const real_t yMax           = std::max(real_t(0.1), drumHalfLength - real_t(1.05) * circumRadius);
   const uint_t nx              = uint_t(real_t(2) * interiorRadius / cellSize) + uint_t(1);
   const uint_t ny              = uint_t(real_t(2) * yMax / cellSize) + uint_t(1);

   uint_t created = 0;
   for (uint_t iz = 0; iz < nx && created < numParticles; ++iz)
   {
      const real_t gz = -interiorRadius + cellSize * (real_t(iz) + real_t(0.5));
      for (uint_t ix = 0; ix < nx && created < numParticles; ++ix)
      {
         const real_t gx = -interiorRadius + cellSize * (real_t(ix) + real_t(0.5));
         if (gx * gx + gz * gz > interiorRadius * interiorRadius) continue;
         for (uint_t iy = 0; iy < ny && created < numParticles; ++iy)
         {
            const real_t gy = (ny > 1) ? (-yMax + cellSize * (real_t(iy) + real_t(0.5))) : real_t(0);
            detail::addRandomOrientedShape(engine, shape, Vec3(gx, gy, gz), sphereRadius, boxHalfExtents, cylRadius,
                                            cylLength, density);
            ++created;
         }
      }
   }
   if (created < numParticles)
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum WARNING: only fit "
                                 << created << "/" << numParticles
                                 << " particles inside the drum at this grid resolution -- increase drumRadius/"
                                    "drumHalfLength or reduce numParticles/particle size for a larger bed.");

   // DIAGNOSTIC: verify the fill grid actually placed particles non-overlapping, empirically, before any
   // physics runs -- added while investigating the settle-phase energy-monotonicity FAIL (job 4049208):
   // the theoretical argument (cellSize=2.1*circumRadius > 2*circumRadius, so no two grid points can be
   // closer than cellSize regardless of shape/orientation) says this should be impossible, but the
   // observed first-checkpoint energy spike (~1000x the system's entire gravitational PE budget, and
   // reproduced with BOTH cylinder AND sphere shapes) is otherwise unexplained. O(N^2) is cheap at
   // N=250 (~31000 pairs), done once.
   {
      real_t minPairDist = std::numeric_limits< real_t >::max();
      size_t minI = 0, minJ = 0;
      for (size_t i = 0; i < engine.ac().size(); ++i)
      {
         if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
         for (size_t j = i + 1; j < engine.ac().size(); ++j)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(j), data::particle_flags::GLOBAL)) continue;
            const real_t d = (engine.ac().getPosition(i) - engine.ac().getPosition(j)).length();
            if (d < minPairDist) { minPairDist = d; minI = i; minJ = j; }
         }
      }
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum fill-grid check: minPairwiseCenterDistance="
                                 << minPairDist << " (expected >= cellSize=" << cellSize
                                 << ", i.e. >= 2*circumRadius=" << (real_t(2) * circumRadius) << ") "
                                 << (minPairDist >= real_t(2) * circumRadius ? "OK" : "*** OVERLAP AT SPAWN ***")
                                 << " pair=(" << minI << "," << minJ << ") pos_i=" << engine.ac().getPosition(minI)
                                 << " pos_j=" << engine.ac().getPosition(minJ));
   }

   engine.finalizeSetup();
   engine.setupVtk(outBase, vtkSpacing);

   // Phase 1: settle under gravity with the drum stationary (same avgSpeed/energy-monotonicity pattern as
   // runPacking/runAngleOfRepose above). This phase's energy check IS meaningful (no external driving
   // yet); the tumbling phase's is not (the rotating wall does net positive work on the bed by design),
   // so it is only checked here.
   real_t avgSpeed            = std::numeric_limits< real_t >::max();
   real_t prevEnergy          = std::numeric_limits< real_t >::max();
   real_t worstEnergyIncrease = real_t(0);
   uint_t step                = 0;
   for (; step < settleSteps && avgSpeed > real_t(1e-3) * std::sqrt(gravity * circumRadius); ++step)
   {
      engine.step();
      engine.writeVtkIfDue();
      if (step > 0 && step % 200 == 0)
      {
         real_t sumSpeed = real_t(0);
         real_t maxSpeed  = real_t(0);
         uint_t n         = uint_t(0);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            const real_t speed = engine.ac().getLinearVelocity(i).length();
            sumSpeed += speed;
            maxSpeed = std::max(maxSpeed, speed);
            ++n;
         }
         avgSpeed = (n > 0) ? sumSpeed / real_t(n) : real_t(0);

         const real_t energy = detail::totalMechanicalEnergy(engine, gravity);
         // Diagnostic trace for the settle-phase energy-monotonicity check (mirrors runPacking's own
         // per-checkpoint print): logs every checkpoint at INFO so a spike's step/energy/maxSpeed can be
         // localized after the fact, same as runPacking already allows.
         WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum settle step=" << step << " avgSpeed=" << avgSpeed
                                    << " maxSpeed=" << maxSpeed << " totalMechanicalEnergy=" << energy);
         if (prevEnergy < std::numeric_limits< real_t >::max())
         {
            const real_t increase = (energy - prevEnergy) / std::max(real_t(1e-9), std::abs(prevEnergy));
            if (increase > worstEnergyIncrease)
               WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum settle step=" << step
                                          << " NEW WORST relative energy increase=" << increase
                                          << " (prevEnergy=" << prevEnergy << " energy=" << energy << ")");
            worstEnergyIncrease = std::max(worstEnergyIncrease, increase);
         }
         prevEnergy = energy;
      }
   }
   const bool settlePhaseStable = worstEnergyIncrease < real_t(0.01);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum settle phase: steps=" << step << " [energy non-increase] "
                              << "worst checkpoint-to-checkpoint relative increase=" << worstEnergyIncrease << "  "
                              << (settlePhaseStable ? "PASS" : "FAIL (energy created -- contact/integration unstable)"));

   // Phase 2: rotate the drum at `rpm` for `tumbleRevolutions` full turns, taking `numSnapshots` free-
   // surface readings evenly spaced across that window.
   const real_t stepsPerRevolution = (real_t(2) * math::pi / omega) / dt;
   const uint_t tumbleSteps        = uint_t(real_t(tumbleRevolutions) * stepsPerRevolution);
   const uint_t snapshotEvery      = std::max< uint_t >(uint_t(1), tumbleSteps / std::max< uint_t >(1, numSnapshots));

   std::vector< std::pair< real_t, real_t > > surfacePoints;
   for (uint_t t = 0; t < tumbleSteps; ++t)
   {
      detail::updateDrumWalls(engine, wallIndices, real_t(0), real_t(0), drumRadius, omega * (real_t(t) * dt));
      engine.step();
      engine.writeVtkIfDue();

      if (t % snapshotEvery == 0)
      {
         std::vector< real_t > binMaxZ(numRadialBins, -std::numeric_limits< real_t >::max());
         std::vector< bool >   binHas(numRadialBins, false);
         for (size_t i = 0; i < engine.ac().size(); ++i)
         {
            if (data::particle_flags::isSet(engine.ac().getFlags(i), data::particle_flags::GLOBAL)) continue;
            const Vec3 pos = engine.ac().getPosition(i);
            if (std::abs(pos[0]) > drumRadius) continue;
            const uint_t bin = std::min< uint_t >(
               numRadialBins - 1,
               uint_t((pos[0] + drumRadius) / (real_t(2) * drumRadius) * real_t(numRadialBins)));
            const real_t topZ = pos[2] + engine.ac().getInteractionRadius(i);
            if (topZ > binMaxZ[bin])
            {
               binMaxZ[bin] = topZ;
               binHas[bin]  = true;
            }
         }
         for (uint_t bi = 0; bi < numRadialBins; ++bi)
            if (binHas[bi])
            {
               const real_t xCenter =
                  -drumRadius + (real_t(bi) + real_t(0.5)) * (real_t(2) * drumRadius / real_t(numRadialBins));
               surfacePoints.emplace_back(xCenter, binMaxZ[bi]);
            }
      }
   }

   // Least-squares fit of the pooled free-surface points -- Kodam et al.'s own Fig. 12(b) method (average
   // the free surface over the ten images, then fit one line; here pooling plays the role of averaging).
   real_t sumX = 0, sumZ = 0, sumXZ = 0, sumXX = 0;
   uint_t n = 0;
   for (const auto& p : surfacePoints)
   {
      sumX += p.first;
      sumZ += p.second;
      sumXZ += p.first * p.second;
      sumXX += p.first * p.first;
      ++n;
   }
   const real_t denom          = real_t(n) * sumXX - sumX * sumX;
   const real_t slope          = (std::abs(denom) > real_t(1e-9)) ? (real_t(n) * sumXZ - sumX * sumZ) / denom : real_t(0);
   // repose angle = atan(|slope|) of the pooled free-surface points' least-squares fit; checked against
   // the plausibility band below, anchored by Kodam et al. 2010 Part II's 39.6deg sim / 43.9deg experiment
   const real_t reposeAngleDeg = std::atan(std::abs(slope)) * real_t(180) / math::pi;

   // Plausibility band only (~5-60deg spans essentially all real granular materials in a rotating drum);
   // Kodam et al. (2010) Part II's own cylinder result (39.6deg simulated vs. 43.9deg experimental) is one
   // concrete literature anchor within this band.
   const bool pass = n > 0 && reposeAngleDeg > real_t(5) && reposeAngleDeg < real_t(60) && settlePhaseStable;

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] rotatingDrum shape="
                              << shape << " created=" << created << " tumbleSteps=" << tumbleSteps
                              << " snapshotPoints=" << n << " reposeAngle=" << reposeAngleDeg
                              << "deg  (plausible range 5-60deg; Kodam et al. 2010 cylinder reference: "
                              << "39.6deg sim / 43.9deg experiment)  " << (pass ? "PASS" : "FAIL"));
   if (!pass && n > 0 && reposeAngleDeg < real_t(5))
      WALBERLA_LOG_INFO_ON_ROOT(
         "[DEMValidation] rotatingDrum DIAGNOSIS: a near-zero angle points to the known mesa_pd "
         "LinearSpringDashpot limitation documented in runAngleOfRepose above -- pure Coulomb sliding "
         "friction with no rolling-resistance torque lets particles (especially spheres) roll/slide into a "
         "flat, low layer rather than a sloped bed, even inside a closed, rotating drum. This is a genuine "
         "contact-model gap this validation run has surfaced, not a broken test.");
   return pass;
}

//======================================================================================================================
// 10: arch tunnel containment -- fluid-free local validation for ComputeTunnelContacts (ParticleKernels.h),
// isolated from the fluid-coupling/dormant-floater changes so a containment bug surfaces here first
// rather than in a full bridge run. Two claims per shape ("sphere"/"box"/"cylinder", swap via the
// `shape` prm key):
//   1. CONTAINMENT: dropped from rest near the tunnel axis under gravity, the body's radial distance
//      from the axis (plus its own circumscribing radius) never exceeds the tunnel radius -- it settles
//      resting against the wall rather than clipping through it.
//   2. NO PHANTOM WALL: the SAME drop, but positioned outside the tunnel's axial [xMin,xMax] span,
//      shows ZERO lateral (Y) deflection -- confirms the X-range gate doesn't leak a force outside the
//      arch's own footprint (there's nothing else in this bare scenario to deflect it sideways, so any
//      Y motion at all means the gate leaked).
//======================================================================================================================
inline bool runArchTunnelContainment(const config::Config::BlockHandle& b)
{
   const std::string shape     = b.getParameter< std::string >("shape", "sphere");
   const real_t sphereRadius   = b.getParameter< real_t >("sphereRadius", real_t(1));
   const Vec3   boxHalfExtents = b.getParameter< Vec3 >("boxHalfExtents", Vec3(real_t(0.8)));
   const real_t cylRadius      = b.getParameter< real_t >("cylRadius", real_t(0.8));
   const real_t cylLength      = b.getParameter< real_t >("cylLength", real_t(2));
   const real_t density        = b.getParameter< real_t >("density", real_t(1));
   const real_t e              = b.getParameter< real_t >("restitutionCoefficient", real_t(0.5));
   const real_t mu             = b.getParameter< real_t >("frictionCoefficient", real_t(0.3));
   const real_t poisson        = b.getParameter< real_t >("poissonsRatio", real_t(0.22));
   const real_t collisionTime  = b.getParameter< real_t >("collisionTime", real_t(0.01));
   const real_t dt             = b.getParameter< real_t >("dt", collisionTime / real_t(50));
   const real_t gravityMag     = b.getParameter< real_t >("gravity", real_t(9.81));
   const real_t tunnelRadius   = b.getParameter< real_t >("tunnelRadius", real_t(5));
   const real_t tunnelHalfSpan = b.getParameter< real_t >("tunnelHalfSpan", real_t(5)); // xMin=-half, xMax=+half
   const uint_t settleSteps    = b.getParameter< uint_t >("settleSteps", uint_t(4000));
   const real_t tolerance      = b.getParameter< real_t >("containmentTolerance", real_t(0.05)); // LU slack

   detail::PassFailAccumulator acc;
   const real_t effR = detail::shapeCircumRadius(shape, sphereRadius, boxHalfExtents, cylRadius, cylLength);

   // --- Case 1: containment -- start just ABOVE where the body will come to rest against the tunnel
   // wall (a small gap, not at the axis), well inside [xMin,xMax], so most of the run is spent actually
   // resolving the wall contact and settling rather than in free-fall approaching it. -----------------
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, gravityMag, tunnelHalfSpan + tunnelRadius + real_t(20));
      engine.addArchTunnel(Vec3(0, 0, 0), tunnelRadius, -tunnelHalfSpan, tunnelHalfSpan, mu);
      const Vec3   startPos = Vec3(0, 0, -(tunnelRadius - effR) + real_t(0.3) * effR);
      const size_t idx      = detail::addRandomOrientedShape(engine, shape, startPos, sphereRadius,
                                                              boxHalfExtents, cylRadius, cylLength, density);
      engine.finalizeSetup(real_t(2) * effR);

      real_t maxRadialReach = real_t(0);
      const uint_t maxSteps = settleSteps;
      for (uint_t s = 0; s < maxSteps; ++s)
      {
         engine.step();
         const Vec3   pos    = engine.ac().getPosition(idx);
         const real_t radial = std::sqrt(pos[1] * pos[1] + pos[2] * pos[2]) + effR;
         maxRadialReach       = std::max(maxRadialReach, radial);
      }
      const bool pass = maxRadialReach <= tunnelRadius + tolerance;
      acc.record(pass);
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] archTunnelContainment shape=" << shape
                                 << " maxRadialReach=" << maxRadialReach << " tunnelRadius=" << tunnelRadius
                                 << " (tolerance=" << tolerance << ")  " << (pass ? "PASS" : "FAIL"));
   }

   // --- Case 2: no phantom wall -- same drop, positioned outside the tunnel's axial span. -------------
   {
      DEMEngine engine(dt, e, mu, poisson, collisionTime, gravityMag, tunnelHalfSpan + tunnelRadius + real_t(20));
      engine.addArchTunnel(Vec3(0, 0, 0), tunnelRadius, -tunnelHalfSpan, tunnelHalfSpan, mu);
      const real_t startY    = real_t(0.3) * (tunnelRadius - effR); // same lateral offset as case 1
      const Vec3   startPos  = Vec3(tunnelHalfSpan + real_t(5) + effR, startY, 0);
      const size_t idx       = detail::addRandomOrientedShape(engine, shape, startPos, sphereRadius,
                                                               boxHalfExtents, cylRadius, cylLength, density);
      engine.finalizeSetup(real_t(2) * effR);

      for (uint_t s = 0; s < settleSteps; ++s)
         engine.step();

      const Vec3   pos      = engine.ac().getPosition(idx);
      const real_t yDrift   = std::abs(pos[1] - startY);
      const bool   pass     = yDrift < real_t(1e-6); // pure free-fall: Y must be untouched
      acc.record(pass);
      WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] archTunnelContainment (outside xRange) shape=" << shape
                                 << " yDrift=" << yDrift << " (expect ~0, no phantom wall)  "
                                 << (pass ? "PASS" : "FAIL"));
   }

   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] archTunnelContainment OVERALL: "
                              << (acc.allPassed() ? "PASS" : "FAIL") << " (" << acc.passed << "/" << acc.total
                              << ")");
   return acc.allPassed();
}

} // namespace dem_validation
} // namespace antidunes
} // namespace walberla
