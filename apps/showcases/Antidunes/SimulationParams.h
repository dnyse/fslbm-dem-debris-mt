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
//! \file SimulationParams.h
//! \author Behzad Safaei <behzad.safaei@fau.de>
//
// Plain-data parameter structs. Nothing in here owns any waLBerla resources.
//======================================================================================================================

#pragma once

#include "core/config/Config.h"
#include "core/math/Vector3.h"

#include "Utility.h"


namespace walberla {
namespace antidunes {

struct BlockForestParams {
   Vector3< uint_t > cellsPerBlock;
   Vector3< bool >   periodicity;
   bool              loadSnapshot;
   bool              storeSnapshot;
   uint_t            snapshotFrequency;
   std::string       snapshotBaseFolder;

   BlockForestParams(const config::Config::BlockHandle& b)
      : loadSnapshot       { b.getParameter< bool >               ("loadSnapshot") }
      , storeSnapshot      { b.getParameter< bool >               ("storeSnapshot") }
      , snapshotFrequency  { b.getParameter< uint_t >             ("snapshotFrequency") }
      , snapshotBaseFolder { b.getParameter< std::string >        ("snapshotBaseFolder") }
   {
      cellsPerBlock =        b.getParameter< Vector3< uint_t > >  ("cellsPerBlock");
      periodicity =          b.getParameter< Vector3< bool > >    ("periodicity");
   }
};

struct DomainParams {
   Vector3< uint_t > domainSize;
   uint_t            wavePeriods;
   real_t            liquidHeightFactor;
   real_t            floorHeightFactor;
   real_t            initialAmplitude;

   DomainParams(const config::Config::BlockHandle& b)
      : wavePeriods         { b.getParameter< uint_t >            ("wavePeriods") }
      , liquidHeightFactor  { b.getParameter< real_t >            ("liquidHeightFactor") }
      , floorHeightFactor   { b.getParameter< real_t >            ("floorHeightFactor") }
      , initialAmplitude    { b.getParameter< real_t >            ("initialAmplitude") }
   {
      domainSize =            b.getParameter< Vector3< uint_t > > ("domainSize");
   }
};

struct PIDParams {
   real_t targetMeanVelocityMagnitude;
   real_t proportionalGain;
   real_t derivativeGain;
   real_t integralGain;
   real_t maxRamp;
   real_t minActuatingVariable;
   real_t maxActuatingVariable;

   PIDParams(const config::Config::BlockHandle& b)
      : targetMeanVelocityMagnitude { b.getParameter< real_t >("targetMeanVelocityMagnitude") }
      , proportionalGain            { b.getParameter< real_t >("proportionalGain") }
      , derivativeGain              { b.getParameter< real_t >("derivativeGain") }
      , integralGain                { b.getParameter< real_t >("integralGain") }
      , maxRamp                     { b.getParameter< real_t >("maxRamp") }
      , minActuatingVariable        { b.getParameter< real_t >("minActuatingVariable") }
      , maxActuatingVariable        { b.getParameter< real_t >("maxActuatingVariable") }
   {}
};

struct PhysicsParams {
   bool   enableWetting;
   uint_t timesteps;
   real_t Re;
   real_t Fr;
   real_t We;

   PhysicsParams(const config::Config::BlockHandle& b)
      : enableWetting { b.getParameter< bool >  ("enableWetting") }
      , timesteps     { b.getParameter< uint_t >("timesteps") }
      , Re            { b.getParameter< real_t >("Re") }
      , Fr            { b.getParameter< real_t >("Fr") }
      , We            { b.getParameter< real_t >("We") }
   {}
};

struct ParticleParams {
   std::string inFileName;
   uint_t      bedCopiesInX;
   uint_t      bedCopiesInY;
   real_t      densityRatio;
   real_t      fixingHeightFactor;
   real_t      frictionCoefficient;
   real_t      restitutionCoefficient;
   uint_t      numSubCycles;
   bool        useLubricationCorrection;
   bool        useNoSlipParticles;
   bool        useParticleBed;
   real_t      floaterRingBand;   // sphere V_sub plane-fit: keep interface cells with |dist-r|<band.
                                  // 1.5 = contact-ring restriction (default). A large value (e.g. 1e6)
                                  // reproduces the old full-box sampling, for clean A/B comparison.

   ParticleParams(const config::Config::BlockHandle& b)
      : inFileName               { b.getParameter< std::string >("inFileName") }
      , bedCopiesInX             { b.getParameter< uint_t >     ("bedCopiesInX") }
      , bedCopiesInY             { b.getParameter< uint_t >     ("bedCopiesInY") }
      , densityRatio             { b.getParameter< real_t >     ("densityRatio") }
      , fixingHeightFactor       { b.getParameter< real_t >     ("fixingHeightFactor") }
      , frictionCoefficient      { b.getParameter< real_t >     ("frictionCoefficient") }
      , restitutionCoefficient   { b.getParameter< real_t >     ("restitutionCoefficient") }
      , numSubCycles             { b.getParameter< uint_t >     ("numSubCycles") }
      , useLubricationCorrection { b.getParameter< bool >       ("useLubricationCorrection") }
      , useNoSlipParticles       { b.getParameter< bool >       ("useNoSlipParticles") }
      , useParticleBed           { b.getParameter< bool >       ("useParticleBed", true) }
      , floaterRingBand          { b.getParameter< real_t >     ("floaterRingBand", real_t(1.5)) }
   {}
};

struct PillarParams {
   Vector3< real_t > center;
   Vector3< real_t > halfExtents;

   PillarParams(const config::Config::BlockHandle& b)
   {
      center      = b.getParameter< Vector3< real_t > >("center");
      halfExtents = b.getParameter< Vector3< real_t > >("halfExtents");
   }
};

struct CylinderPillarParams {
   Vector3< real_t > center;
   real_t            radius;

   CylinderPillarParams(const config::Config::BlockHandle& b)
   {
      center = b.getParameter< Vector3< real_t > >("center");
      radius = b.getParameter< real_t >           ("radius");
   }
};

struct FloatingSphereParams {
   Vector3< real_t > position;
   real_t            radius;
   real_t            densityRatio;
   // Number of LBM timesteps the sphere is held fixed before being released to float freely.
   // Use this when starting with the sphere at or above the free surface: the FSLBM needs time
   // to establish a stable fill-level distribution around the contact line before the sphere moves.
   // 0 = release immediately (default, correct for fully submerged start).
   uint_t            settlingSteps;

   FloatingSphereParams(const config::Config::BlockHandle& b)
   {
      position      = b.getParameter< Vector3< real_t > >("position");
      radius        = b.getParameter< real_t >           ("radius");
      densityRatio  = b.getParameter< real_t >           ("densityRatio");
      settlingSteps = b.getParameter< uint_t >           ("settlingSteps", uint_t(0));
   }
};

inline std::vector< PillarParams > parsePillars(const shared_ptr< config::Config >& cfg)
{
   std::vector< PillarParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("Pillars", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks pillarBlocks;
   containers[0].getBlocks("Pillar", pillarBlocks, 0, 100);
   for (const auto& b : pillarBlocks)
      result.emplace_back(b);
   return result;
}

inline std::vector< CylinderPillarParams > parseCylinderPillars(const shared_ptr< config::Config >& cfg)
{
   std::vector< CylinderPillarParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("CylinderPillars", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks pillarBlocks;
   containers[0].getBlocks("CylinderPillar", pillarBlocks, 0, 100);
   for (const auto& b : pillarBlocks)
      result.emplace_back(b);
   return result;
}

inline std::vector< FloatingSphereParams > parseFloatingSpheres(const shared_ptr< config::Config >& cfg)
{
   std::vector< FloatingSphereParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("FloatingSpheres", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks sphereBlocks;
   containers[0].getBlocks("FloatingSphere", sphereBlocks, 0, 100);
   for (const auto& b : sphereBlocks)
      result.emplace_back(b);
   return result;
}

struct FloatingBoxParams {
   Vector3< real_t > position;
   Vector3< real_t > halfExtents;   // half-size in each dimension; edgeLength = 2 * halfExtents
   real_t            densityRatio;
   uint_t            settlingSteps;      // same semantics as FloatingSphereParams::settlingSteps
   real_t            initialTiltAngleDeg; // initial rotation about the Y-axis (degrees); 0 = upright

   FloatingBoxParams(const config::Config::BlockHandle& b)
   {
      position             = b.getParameter< Vector3< real_t > >("position");
      halfExtents          = b.getParameter< Vector3< real_t > >("halfExtents");
      densityRatio         = b.getParameter< real_t >           ("densityRatio");
      settlingSteps        = b.getParameter< uint_t >           ("settlingSteps", uint_t(0));
      initialTiltAngleDeg  = b.getParameter< real_t >           ("initialTiltAngleDeg", real_t(0));
   }
};

inline std::vector< FloatingBoxParams > parseFloatingBoxes(const shared_ptr< config::Config >& cfg)
{
   std::vector< FloatingBoxParams > result;
   config::Config::Blocks containers;
   cfg->getBlocks("FloatingBoxes", containers, 0, 1);
   if (containers.empty()) return result;
   config::Config::Blocks boxBlocks;
   containers[0].getBlocks("FloatingBox", boxBlocks, 0, 100);
   for (const auto& b : boxBlocks)
      result.emplace_back(b);
   return result;
}

struct ModelParams {
   std::string pdfReconstructionModel;
   std::string pdfRefillingModel;
   std::string excessMassDistributionModel;
   std::string curvatureModel;
   bool        useSimpleMassExchange;
   real_t      cellConversionThreshold;
   real_t      cellConversionForceThreshold;
   // true  → Bogner & Rüde (2013) per-cell gas/liquid/interface state selection for FormerMO cells.
   // false → legacy uniform reconstruction: every uncovered cell becomes an interface cell
   //         (no fill=1 over-promotion band; the all-interface path that never drifted/NaN'd).
   bool        useBognerReconstruction;

   ModelParams(const config::Config::BlockHandle& b)
      : pdfReconstructionModel      { b.getParameter< std::string >("pdfReconstructionModel") }
      , pdfRefillingModel           { b.getParameter< std::string >("pdfRefillingModel") }
      , excessMassDistributionModel { b.getParameter< std::string >("excessMassDistributionModel") }
      , curvatureModel              { b.getParameter< std::string >("curvatureModel") }
      , useSimpleMassExchange       { b.getParameter< bool >       ("useSimpleMassExchange") }
      , cellConversionThreshold     { b.getParameter< real_t >     ("cellConversionThreshold") }
      , cellConversionForceThreshold{ b.getParameter< real_t >     ("cellConversionForceThreshold") }
      , useBognerReconstruction     { b.getParameter< bool >       ("useBognerReconstruction", true) }
   {}
};

struct EvaluationParams {
   uint_t      performanceLogFrequency;
   uint_t      evaluationFrequency;
   std::string baseFolderName;

   EvaluationParams(const config::Config::BlockHandle& b)
      : performanceLogFrequency { b.getParameter< uint_t >     ("performanceLogFrequency") }
      , evaluationFrequency     { b.getParameter< uint_t >     ("evaluationFrequency") }
      , baseFolderName          { b.getParameter< std::string >("baseFolderName") }
   {}
};

struct VTKParams {
   uint_t            writeFrequency;
   uint_t            ghostLayers;
   std::string       baseFolder;
   real_t            samplingResolution;

   VTKParams(const config::Config::BlockHandle& b)
      : writeFrequency    { b.getOneBlock("fluid_field").getParameter< uint_t >     ("writeFrequency") }
      , ghostLayers       { b.getOneBlock("fluid_field").getParameter< uint_t >     ("ghostLayers") }
      , baseFolder        { b.getOneBlock("fluid_field").getParameter< std::string >("baseFolder") }
      , samplingResolution{ b.getOneBlock("fluid_field").getParameter< real_t >     ("samplingResolution") }
   {}
};

struct SimulationParams;

// Parameters that are computed once and for all (derived from config params)
//-----------------------------------------------------------------------
struct ComputedParams {
   real_t avgParticleDiameter;
   real_t particleScalingFactor;
   real_t particleKappa;
   real_t particleCollisionTime;
   real_t particleTimestepSize;
   real_t liquidHeight;
   real_t floorHeight;
   real_t absoluteLiquidHeight;
   real_t particleFixingHeight;
   real_t viscosity;
   real_t relaxationRate;
   real_t gravity;
   real_t surfaceTension;
   real_t dx_SI;
   real_t dt_SI;
   real_t sinusAmplitude;
   real_t sinusOffset;
   real_t sinusWavelength;

   ComputedParams(
      const ParticleParams& particle, 
      const DomainParams& domain, 
      const PIDParams& pid, 
      const PhysicsParams& physics)
   {
      getAvgDiameterScalingFactor(
         particle.inFileName, domain.domainSize, 
         particle.bedCopiesInX, particle.bedCopiesInY,
         avgParticleDiameter, particleScalingFactor);

      const real_t particlePoissonsRatio = real_c(0.22);
      particleKappa = real_c(2) * (real_c(1) - particlePoissonsRatio) / (real_c(2) - particlePoissonsRatio);
      const real_t particleCollisionTimeNonDim = real_c(4);
      particleCollisionTime = particleCollisionTimeNonDim * avgParticleDiameter;
      particleTimestepSize = real_c(1) / real_c(particle.numSubCycles);
      
      liquidHeight         = avgParticleDiameter * domain.liquidHeightFactor;
      floorHeight          = avgParticleDiameter * domain.floorHeightFactor;
      WALBERLA_CHECK_FLOAT_UNEQUAL(liquidHeight, 0.0)
      WALBERLA_CHECK_FLOAT_UNEQUAL(floorHeight, 0.0)
      absoluteLiquidHeight = liquidHeight + floorHeight;
      particleFixingHeight = particle.fixingHeightFactor * avgParticleDiameter;

      const real_t U  = pid.targetMeanVelocityMagnitude;
      viscosity      = U * liquidHeight / physics.Re;                            // 0.000214026
      relaxationRate = real_c(1) / (real_c(3) * viscosity + real_c(0.5));        // 1.99743
      gravity        = (U / physics.Fr) * (U / physics.Fr) / liquidHeight;       // 7.02618e-06
      surfaceTension = real_c(1) * U * U * liquidHeight / physics.We;            // 0.000849592

      const real_t viscosity_SI = real_c(1.0016e-6);  // kinemtic viscosity of water at 20°C at 1 bar
      dx_SI = real_c(1) / particleScalingFactor;
      dt_SI = viscosity / viscosity_SI * dx_SI * dx_SI;

      sinusAmplitude  = real_c(0.5) * domain.initialAmplitude;
      sinusOffset     = floorHeight;
      sinusWavelength = real_c(domain.domainSize[0]) / real_c(domain.wavePeriods);
   }

   void logOnRoot() const
   {
      WALBERLA_LOG_INFO_ON_ROOT("Computed parameters --------------------")
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(avgParticleDiameter);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleScalingFactor);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleKappa);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleCollisionTime);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleTimestepSize);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(liquidHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(floorHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(absoluteLiquidHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleFixingHeight);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationRate);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(gravity);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(surfaceTension);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dx_SI);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dt_SI);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(sinusAmplitude);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(sinusOffset);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(sinusWavelength);
      WALBERLA_LOG_INFO_ON_ROOT("------------------------------------------")
   }

   void writeToFile(const std::string& baseFolderName) const
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file(baseFolderName + "/computedParams.txt");
         file << "avgParticleDiameter:   " << avgParticleDiameter    << "\n";
         file << "particleScalingFactor: " << particleScalingFactor  << "\n";
         file << "particleKappa:         " << particleKappa          << "\n";
         file << "particleCollisionTime: " << particleCollisionTime  << "\n";
         file << "particleTimestepSize:  " << particleTimestepSize   << "\n";
         file << "liquidHeight:          " << liquidHeight           << "\n";
         file << "floorHeight:           " << floorHeight            << "\n";
         file << "absoluteLiquidHeight:  " << absoluteLiquidHeight   << "\n";
         file << "particleFixingHeight:  " << particleFixingHeight   << "\n";
         file << "viscosity:             " << viscosity              << "\n";
         file << "relaxationRate:        " << relaxationRate         << "\n";
         file << "gravity:               " << gravity                << "\n";
         file << "surfaceTension:        " << surfaceTension         << "\n";
         file << "dx_SI:                 " << dx_SI                  << "\n";
         file << "dt_SI:                 " << dt_SI                  << "\n";
         file << "sinusAmplitude:        " << sinusAmplitude         << "\n";
         file << "sinusOffset:           " << sinusOffset            << "\n";
         file << "sinusWavelength:       " << sinusWavelength        << "\n";
         file.close();
      }
   }
};

// Container for all simulation parameters 
//-----------------------------------------------------------------------
struct SimulationParams {
   // config ----------
   shared_ptr< config::Config > cfgPtr;
   
   // config params ----------
   BlockForestParams blockForest;
   DomainParams      domain;
   PIDParams         pid;
   PhysicsParams     physics;
   ParticleParams                       particle;
   std::vector< PillarParams >          pillars;
   std::vector< CylinderPillarParams >  cylinderPillars;
   std::vector< FloatingSphereParams >  floatingSpheres;
   std::vector< FloatingBoxParams >     floatingBoxes;
   ModelParams       model;
   EvaluationParams  evaluation;
   VTKParams         vtk;

   // computed params ----------
   ComputedParams    computed;

   SimulationParams(const shared_ptr< config::Config >& cfg)
      : cfgPtr       { cfg }
      , blockForest  { cfg->getOneBlock("BlockForestParameters") }
      , domain       { cfg->getOneBlock("DomainParameters") }
      , pid          { cfg->getOneBlock("PIDParameters") }
      , physics      { cfg->getOneBlock("PhysicsParameters") }
      , particle     { cfg->getOneBlock("ParticleParameters") }
      , pillars          { parsePillars(cfg) }
      , cylinderPillars  { parseCylinderPillars(cfg) }
      , floatingSpheres  { parseFloatingSpheres(cfg) }
      , floatingBoxes    { parseFloatingBoxes(cfg) }
      , model        { cfg->getOneBlock("ModelParameters") }
      , evaluation   { cfg->getOneBlock("EvaluationParameters") }
      , vtk          { cfg->getOneBlock("VTK") }
      , computed     { particle, domain, pid, physics }
   {}
};


} // namespace antidunes
} // namespace walberla