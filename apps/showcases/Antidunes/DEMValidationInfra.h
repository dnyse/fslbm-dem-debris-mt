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
//! wall-impact loop, energy/conservation logging helpers), included by DEMValidationWallImpact.h,
//! DEMValidationPairEnergy.h and DEMValidationPackingAndDrum.h.

#pragma once

#include <sstream>

#include "AnalyticalSolutions.h"
#include "Utility.h"

#include "blockforest/Initialization.h"
// ParticleKernels.h -> SimulationObjects.h uses blockforest::SimpleCommunication<> but relies on the caller
// including its defining header (in Antidunes.cpp via BuildTimeloop.h's lbm/free_surface chain, which this
// fluid-free app does not include) -- included explicitly here.
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
//! sub-cycle at a time. Always uses LinkedCells + forEachParticlePairHalf (as BedGeneration.cpp) regardless of
//! particle count, so one code path serves the 2-body wall-impact tests and the packing tests alike.
class DEMEngine
{
public:
   // useAnalyticalContactDetection: false (default) routes Box x Box / Cylinder x Cylinder through the generic
   // GJK/EPA fallback; true routes them through the production analytical dispatch (boxBoxSAT / Kodam+Guo).
   // GJK/EPA and the closed forms are different code with different failure modes, so a defect specific to one
   // is only reproduced by running that one.
   DEMEngine(real_t dt, real_t restitutionCoefficient, real_t frictionCoefficient, real_t poissonsRatio,
             real_t collisionTime, real_t gravityMagnitude, real_t domainHalfExtent = real_t(500),
             bool useAnalyticalContactDetection = false, bool useRotationalEffectiveMass = false)
      : dt_(dt)
      , gravityZ_(-gravityMagnitude)
      , restitutionCoefficient_(restitutionCoefficient)
      , collisionTime_(collisionTime)
      // kappa = kT/kN = 2(1-nu)/(2-nu), Mindlin (1949)/Thornton's no-slip single-contact stiffness ratio. The
      // primary Thornton reference is not verified against a source; read it as "Thornton-type no-slip model".
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
         /*useContactFilterOwnerTiebreak=*/false, /*useTightLayer1=*/false, useRotationalEffectiveMass);
   }

   //! direct access to the production contact functor (detection-mode switches)
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

   //! Re-poses an existing wall (from addWall) in place. HalfSpace's normal is immutable (no setter in
   //! mesa_pd/data/shape/HalfSpace.h), so its baseShape is replaced by a new HalfSpace(newNormal) alongside
   //! setPosition. Animates the polygon-of-planes drum wall (detail::addDrumWalls/updateDrumWalls); cheap enough
   //! every step for the handful of GLOBAL, non-communicating wall particles.
   void setWallPose(size_t idx, const Vec3& pos, const Vec3& normal)
   {
      ps_->setPosition(idx, pos);
      ps_->setBaseShape(idx, std::make_shared< data::HalfSpace >(normal));
   }

   //! Must be called once after all particles are created and before the first step(). `knownMaxDiameter`,
   //! when nonzero, overrides the scan of existing particles below -- required for scenarios (e.g.
   //! runPackingSingleShape) that call finalizeSetup() before any DYNAMIC particles exist (only GLOBAL walls),
   //! where the scan would size the linked-cell grid at ~0 and every later particle would trip
   //! InsertParticleIntoLinkedCells' "interaction radius too large for this cell size" assertion.
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

   //! One DEM sub-cycle: pre-force integration, contact resolution (production ComputeContacts, as in
   //! doDEMSubCycles), gravity, post-force integration. No lubrication, hydrodynamic or buoyancy force.
   void step()
   {
      ps_->forEachParticle(false, kernel::SelectLocal(), *ac_, preVV_, *ac_);
      syncNN_(*ps_, forest_, domain_);

      data::LinkedCells linkedCells(domain_->getUnionOfLocalAABBs().getExtended(linkedCellWidth_),
                                     linkedCellWidth_);
      // LinkedCells's constructor leaves infiniteParticles_ at its default-member-initializer value (0, via
      // std::atomic<int>'s C++20 value-init) instead of the -1 sentinel its traversal expects; only clear() sets
      // it (BedGeneration.cpp calls clear() too). Without it, the first inserted INFINITE particle's next-pointer
      // chains to particle index 0, so a pair traversal reaching the end of the infinite-particle list also
      // revisits index 0 -- a self-pairing crash when index 0 is the sole finite body (the wall-impact scenarios).
      linkedCells.clear();
      kernel::InsertParticleIntoLinkedCells ipilc;
      ps_->forEachParticle(false, kernel::SelectAll(), *ac_, ipilc, *ac_, linkedCells);
      linkedCells.forEachParticlePairHalf(
         false, kernel::ExcludeInfiniteInfinite(), *ac_,
         [this](const size_t idx1, const size_t idx2, auto& ac) { (*computeContacts_)(idx1, idx2, ac); }, *ac_);

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

   //! Query-only contact check (no force applied) with the production contact-detection functor
   //! (SolidCylinderContactDetection via DoubleCast): the "has this pair started/stopped touching" trigger of the
   //! wall-impact / pair-restitution scenarios.
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

   // edgeLength/rotation: the same two fields, in the same convention, as BuildTimeloop.h's box floater VTK output
   // ("AddOutputFloatingBoxes"), so a box renders as its ORIENTED solid in ParaView (Filters > Glyph, Glyph Type =
   // Box, Scale Array = edgeLength, Scale Mode = Scale by Components; then Filters > Tensor Glyph using
   // "rotation", Extract Eigenvalues OFF). Vec3(0) for non-box shapes (a zero-size box glyph does not render;
   // "radius"-based glyphing still works for those).
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
   // Transposed body->world rotation matrix: ParaView's column-major Tensor Glyph needs the transpose to show the
   // orientation, not its inverse (see BuildTimeloop.h). Visualisation only; kineticEnergy() and the production
   // buoyancy/torque math use getRotation().getMatrix() untransposed.
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

//! KE of a rigid body given explicit state, for the ANALYTICAL model's predicted post-impact energy
//! (DEMEngine::kineticEnergy reads a live simulated particle instead).
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

// velocityScale/angularVelocityScale (default 0: pure relative error with the 1e-8 floor): optional mixed
// absolute+relative tolerance, the standard combined form (atol + rtol*|b|) rewritten as a floored-denominator
// relative error, relErr = resid / (scale + |analytic|). The pure relative form is unusable where the analytical
// reference crosses zero: at a tilt where the true rebound speed is ~0.01 of the impact speed, an ordinary small
// absolute residual becomes a 10-200x relative error (see the thesis's wall-impact section). `scale` must be a
// fixed physical reference for the quantity (e.g. the impact speed), not the per-tilt analytic value, or the
// tolerance would be circular.
inline void logCaseResult(const std::string& tag, const Vec3& simVel, const Vec3& analyticVel,
                           const Vec3& simAngVel, const Vec3& analyticAngVel, real_t tolerance,
                           PassFailAccumulator& acc, real_t velocityScale = real_t(0),
                           real_t angularVelocityScale = real_t(0))
{
   const real_t velErr = (simVel - analyticVel).length() /
                          std::max(real_t(1e-8), velocityScale + analyticVel.length());
   const real_t angErr = (simAngVel - analyticAngVel).length() /
                          std::max(real_t(1e-8),
                                   angularVelocityScale + std::max(analyticAngVel.length(), simAngVel.length()));
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
//! or friction coefficients. Checked in EVERY collision scenario (wall-impact and pairwise) as a model-agnostic
//! sanity bound alongside each test's quantitative comparison; a violation would mean the contact law itself
//! injects energy, a more basic problem than a friction-model calibration mismatch. tolerance guards only against
//! integration noise (sub-cycled spring-dashpot contacts are not exactly energy-monotonic step-to-step).
inline void logConservation(const std::string& tag, real_t keBefore, real_t keAfter, PassFailAccumulator& acc,
                             real_t tolerance = real_t(1e-3))
{
   const bool pass = keAfter <= keBefore * (real_t(1) + tolerance) + real_t(1e-9);
   acc.record(pass);
   WALBERLA_LOG_INFO_ON_ROOT("[DEMValidation] " << tag << "  [energy non-increase] KE_before=" << keBefore
                              << " KE_after=" << keAfter << " ratio=" << (keAfter / std::max(real_t(1e-12), keBefore))
                              << "  " << (pass ? "PASS" : "FAIL (energy created -- contact law is unstable)"));
}

//! STRONGER companion to logConservation, for scenarios with a computable post-impact energy (the e=1 elastic
//! case, where it must equal keBefore, or wherever an analytical model predicts the post-state): "did not
//! increase" is one-sided and would pass for far too much or too little dissipation, so this checks the ACTUAL
//! dissipated amount against the EXPECTED one.
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
//! precision regardless of restitution/friction. It holds even for a WRONG friction/restitution calibration, so a
//! violation implicates the contact force's Newton's-third-law symmetry, not parameter tuning.
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
//! system with only dissipative internal forces (inelastic/frictional contacts, no energy source), this must be
//! MONOTONICALLY NON-INCREASING; a violation (checked between sampling checkpoints, not step-to-step, since a
//! single sub-cycle is not exactly energy-monotonic) means the contact/integration scheme injects energy.
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
