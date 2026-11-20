//======================================================================================================================
//! \file DEMValidationProductionCollisionReplay.h
//! ARCHIVED 2026-11-20 from DEMValidationPackingAndDrum.h: replay of one production box-box incident through the
//! 9-way box-contact mitigation configuration set (useRoundedBoxCorners, useMultiPointBoxContact, useBoxOverlapCap,
//! useAxisHysteresis, useContactManifold), all archived. Not compiled; uses the pre-archive DEMEngine signature.
//! See ../ARCHIVE_README.md.
//======================================================================================================================
#pragma once

// Direct replay of a real production incident (job 4199585, bridge "widerstacks" box scenario): one
// box independently confirmed "flying" in ParaView (climbed 31 LU above its own spawn height) against
// its contact partner, both seeded with their exact logged position/velocity/angular-velocity/
// orientation at the step their prolonged-contact episode began, run through the same production
// physics (restitutionCoefficient=0.3 "wet/waterlogged wood", frictionCoefficient=0.5,
// useAnalyticalContactDetection=true).
//
// Measured directly on the real production log (an independent measurement, not this replay): total
// system kinetic energy grew by 536% over this window, on a material whose restitution coefficient
// should LOSE roughly 91% of relative energy in a collision. This scenario checks whether an isolated
// two-body replay (no jam, no fluid coupling) reproduces a comparable gain -- if it does, the defect
// is a genuine pairwise contact-response issue, not something that needs the wider jam context.
//
// Fidelity caveat: best-effort, not bit-exact. LinearSpringDashpot's persistent tangential-spring
// history starts empty at this replay's own step 0, whereas the real contact may have engaged a few
// sub-cycles earlier; a close but not bit-identical energy trajectory is the expected outcome.
//
// Loops the same incident through the standard 9-way mitigation config set (identical rows to
// runMultiPointEdgeCases/runFallOntoPileStability, so all three batteries' results line up in one
// table). Default duration is 10000 sub-cycles, not the real log's own 3540-sub-cycle window: the
// baseline's energy peak doesn't land until sub-cycle ~4072, and the elevated plateau afterward is
// itself part of the signature (an exact-window replay alone understates the real gain, 72% vs. 536%,
// because the contact hasn't finished compressing by sub-cycle 3540).
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
   // "useSpheroBoxContact" row removed -- see archive/unused_contact_methods/.
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
