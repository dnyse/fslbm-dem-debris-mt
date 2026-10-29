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
//! \file DEMValidationInfra.h
//! \brief Shared infrastructure for the DEM validation experiments (DEMEngine, PassFailAccumulator,
//! wall-impact loop, energy/conservation logging helpers) -- split out of DEMValidationScenarios.h
//! 2026-10-29 so each scenario family lives in its own file. See that split's sibling files
//! (DEMValidationWallImpact.h, DEMValidationPairEnergy.h, DEMValidationCornerAndMultiPoint.h,
//! DEMValidationPackingAndDrum.h), all of which include this one.

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
   // and closed_form_contacts::boxBoxSAT/boxBoxMultiPointContacts) -- threaded straight through to
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

} // namespace dem_validation
} // namespace antidunes
} // namespace walberla
