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
//! \file ForcingAdjusterTest.cpp
//! Unit test for walberla::antidunes::ForcingAdjuster's mid-run PID switch-off
//! (PIDParameters.switchOffStep / switchOffAveragingWindow): PID control runs until the switch-off
//! step, then acceleration_x is frozen at the MEAN of the PID's recent commanded values and never
//! touched again -- including across a checkpoint resume before and after the switch-off step.
//
//======================================================================================================================

// Forcing.h is not self-contained: it relies on the include order of Antidunes.cpp (free-surface,
// timeloop and boundary-handling headers arrive via BuildTimeloop.h first). Mirror that order here
// (test-side only, no production file touched) instead of depending on an include-order accident.
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include "sqlite/SQLite.h"
#include "timeloop/all.h"

#include "BuildTimeloop.h"

#include <cmath>
#include <map>
#include <vector>

namespace walberla
{
namespace antidunes
{
namespace
{

constexpr uint_t FREQ = 100;

// bit-exact comparisons on reals without tripping -Wfloat-equal
bool same(real_t a, real_t b) { return !(a < b) && !(b < a); }
#define CHECK_SAME(a, b) WALBERLA_CHECK(same((a), (b)))
#define CHECK_DIFFERENT(a, b) WALBERLA_CHECK(!same((a), (b)))

// Synthetic, oscillating "measured" mean velocity so the PID's commanded force genuinely varies from
// update to update (a constant input would make mean == instantaneous and hide a freeze-at-last-value bug).
real_t syntheticVelocity(uint_t step) { return real_c(0.04) + real_c(0.002) * std::sin(real_c(step) * real_c(0.013)); }

shared_ptr< PIDController > makePid()
{
   return make_shared< PIDController >(real_c(0.04), real_c(7e-7), real_c(2e-4), real_c(1e-6), real_c(2e-4),
                                        real_c(1e-4), real_c(0), real_c(1e-3));
}

struct Run
{
   shared_ptr< Vector3< real_t > > accel;
   std::map< uint_t, real_t > afterCall; // absolute step -> acceleration_x right after that call
};

// Drives the adjuster for calls covering absolute steps [beginStep, endStep).
Run drive(uint_t beginStep, uint_t endStep, uint_t switchOffStep, uint_t window, bool enabled, real_t initialAccel)
{
   Run r;
   r.accel                = make_shared< Vector3< real_t > >(initialAccel, real_c(0), real_c(-1e-5));
   auto meanVelocity      = make_shared< Vector3< real_t > >(real_c(0));
   ForcingAdjuster adjuster(meanVelocity, r.accel, makePid(), FREQ, enabled, beginStep, switchOffStep, window);
   for (uint_t step = beginStep; step < endStep; ++step)
   {
      (*meanVelocity)[0] = syntheticVelocity(step);
      adjuster();
      r.afterCall[step] = (*r.accel)[0];
   }
   return r;
}

// Mean of the PID-commanded values (observable as acceleration_x right after an update call) at the given steps.
real_t meanAt(const Run& r, const std::vector< uint_t >& steps)
{
   real_t sum = real_c(0);
   for (uint_t s : steps) sum += r.afterCall.at(s);
   return sum / real_c(steps.size());
}

// Fresh start, PID until step 2000, window 500: the freeze value is the mean of the PID updates at
// steps 1500..1900 (the last update before the switch is at 1900; entries older than the window drop).
void testFreezeAtWindowMean()
{
   const Run r = drive(0, 3000, 2000, 500, true, real_c(7e-7));

   const real_t expected = meanAt(r, {1500, 1600, 1700, 1800, 1900});
   WALBERLA_CHECK_FLOAT_EQUAL(r.afterCall.at(2000), expected);

   // frozen for the rest of the run, bit-identical
   for (uint_t s = 2000; s < 3000; ++s) CHECK_SAME(r.afterCall.at(s), r.afterCall.at(2000));

   // the PID genuinely moved before the switch, and the mean differs from the last commanded value
   // (otherwise this test could not tell "mean" from "last value")
   WALBERLA_CHECK_GREATER(std::fabs(r.afterCall.at(1900) - r.afterCall.at(1500)), real_c(1e-12));
   WALBERLA_CHECK_GREATER(std::fabs(r.afterCall.at(1900) - expected), real_c(1e-12));

   // and before the switch-off step it was still being adjusted
   CHECK_DIFFERENT(r.afterCall.at(1800), r.afterCall.at(1900));
}

// Same setup with the feature off: the PID must keep adjusting past step 2000.
void testFeatureOffKeepsControlling()
{
   const Run r = drive(0, 3000, 0, 500, true, real_c(7e-7));
   CHECK_DIFFERENT(r.afterCall.at(2500), r.afterCall.at(2900));
}

// Window 0 degenerates to the last commanded value.
void testWindowZeroFreezesAtLastValue()
{
   const Run r = drive(0, 2500, 2000, 0, true, real_c(7e-7));
   WALBERLA_CHECK_FLOAT_EQUAL(r.afterCall.at(2000), r.afterCall.at(1900));
}

// Resume BEFORE the switch-off step: absolute steps (not local call counts) decide when to freeze, and
// only PID updates made in this run enter the mean (resume at 1800 -> updates at 1800 and 1900).
void testResumeBeforeSwitchOff()
{
   const Run r = drive(1800, 2600, 2000, 500, true, real_c(7e-7));
   WALBERLA_CHECK_FLOAT_EQUAL(r.afterCall.at(2000), meanAt(r, {1800, 1900}));
   for (uint_t s = 2000; s < 2600; ++s) CHECK_SAME(r.afterCall.at(s), r.afterCall.at(2000));
}

// Resume AFTER the switch-off step (checkpoint already holds the frozen force): no PID history exists,
// and the loaded acceleration must be left exactly as it is.
void testResumeAfterSwitchOffKeepsLoadedValue()
{
   const real_t loaded = real_c(3.3e-7);
   const Run r         = drive(5000, 6000, 2000, 500, true, loaded);
   for (uint_t s = 5000; s < 6000; ++s) CHECK_SAME(r.afterCall.at(s), loaded);
}

// enabled=false takes precedence: nothing is touched, switchOffStep or not.
void testDisabledIgnoresSwitchOff()
{
   const real_t fixed = real_c(7.8e-7);
   const Run r        = drive(0, 3000, 2000, 500, false, fixed);
   for (uint_t s = 0; s < 3000; ++s) CHECK_SAME(r.afterCall.at(s), fixed);
}

} // namespace
} // namespace antidunes
} // namespace walberla

int main(int argc, char** argv)
{
   walberla::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   using namespace walberla::antidunes;
   testFreezeAtWindowMean();
   testFeatureOffKeepsControlling();
   testWindowZeroFreezesAtLastValue();
   testResumeBeforeSwitchOff();
   testResumeAfterSwitchOffKeepsLoadedValue();
   testDisabledIgnoresSwitchOff();
   return EXIT_SUCCESS;
}
