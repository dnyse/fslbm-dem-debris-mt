#!/usr/bin/env python3
"""Global resistance measurement for a clog, per review.md item 5 ("Use the PID acceleration as a
resistance measurement ... needs no literature, just a momentum balance").

With PID control, in a steady state the driving body force balances wall friction plus the drag of
the pillars and the jam. Writing the momentum balance per unit fluid mass, the EXTRA acceleration a
debris run needs over its own no-debris reference run, times the fluid mass it is accelerating,
estimates the extra drag the jam adds:

    resistance = (a_x - a_x_ref) * mean(totalFluidMass)

both a_x and totalFluidMass are already logged every evaluationFrequency steps to eval/antidunes.sqlite
(WriteData.h's writeFluidInfo) -- no new instrumentation needed, only this post-processing step. In
lattice units, rho ~ 1 near the free surface, so totalFluidMass already plays the role of rho*V_liquid.

For a FIXED-FORCING (PID-off) run, a_x is constant by construction and carries no information; the
review's own substitute is the discharge ratio Q/Q_ref (from the flow-rate probe, ClogMeasurement's
optional `flowRateEnabled`), which plays the same "how much has the jam resisted the flow" role
without needing a_x at all.

This turns the "PID masks the clog" problem (the controller raises the driving force exactly to
compensate for the extra resistance, so a downstream observer never sees the flow slow down) into a
measurement instead of a caveat: it is a single global number, not sensitive to probe placement the
way a local elevation/pressure reading is.

Usage:
    # PID-controlled pair (debris run vs. its own no-debris reference, same forcing regime)
    python3 analyze_clog_resistance.py --mode pid \\
        --debris    path/to/debris_run/eval/antidunes.sqlite \\
        --reference path/to/reference_run/eval/antidunes.sqlite \\
        --steady-from 100000

    # Fixed-forcing pair
    python3 analyze_clog_resistance.py --mode fixed \\
        --debris    path/to/debris_run/eval/antidunes.sqlite \\
        --reference path/to/reference_run/eval/antidunes.sqlite \\
        --steady-from 100000
"""
import argparse
import sqlite3
import statistics
from pathlib import Path


def steady_mean(db_path: str, column: str, steady_from: int) -> tuple[float, float, int]:
    con = sqlite3.connect(db_path)
    rows = con.execute(
        f"SELECT {column} FROM runs WHERE timestep >= ? AND {column} IS NOT NULL ORDER BY timestep",
        (steady_from,),
    ).fetchall()
    con.close()
    if not rows:
        raise SystemExit(f"no rows with timestep >= {steady_from} and non-null {column} in {db_path}")
    vals = [r[0] for r in rows]
    return statistics.mean(vals), statistics.pstdev(vals), len(vals)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=["pid", "fixed"], required=True,
                    help="'pid' uses (a_x - a_x_ref) * mean(totalFluidMass); "
                         "'fixed' uses Q/Q_ref (needs clogUpstreamFlowRate, i.e. flowRateEnabled true)")
    ap.add_argument("--debris", required=True, help="path to the debris run's eval/antidunes.sqlite")
    ap.add_argument("--reference", required=True, help="path to the matching no-debris reference run's eval/antidunes.sqlite")
    ap.add_argument("--steady-from", type=int, required=True,
                     help="timestep from which both runs are considered to have reached steady state "
                          "(inclusive); everything before is discarded as start-up transient")
    ap.add_argument("--flow-probe", default="clogUpstreamFlowRate",
                     help="column to use for --mode fixed (default clogUpstreamFlowRate; "
                          "clogDownstreamFlowRate is the other option)")
    args = ap.parse_args()

    for label, path in [("debris", args.debris), ("reference", args.reference)]:
        if not Path(path).exists():
            raise SystemExit(f"{label} sqlite not found: {path}")

    if args.mode == "pid":
        a_debris, a_debris_std, n_debris = steady_mean(args.debris, "acceleration_x", args.steady_from)
        a_ref, a_ref_std, n_ref = steady_mean(args.reference, "acceleration_x", args.steady_from)
        mass_debris, mass_debris_std, _ = steady_mean(args.debris, "totalFluidMass", args.steady_from)

        resistance = (a_debris - a_ref) * mass_debris

        print(f"debris run    : a_x = {a_debris:.6e} +- {a_debris_std:.2e}  (n={n_debris} samples from t>={args.steady_from})")
        print(f"reference run : a_x = {a_ref:.6e} +- {a_ref_std:.2e}  (n={n_ref} samples from t>={args.steady_from})")
        print(f"debris run    : mean totalFluidMass = {mass_debris:.6e} +- {mass_debris_std:.2e}")
        print()
        print(f"extra driving acceleration (a_x - a_x_ref) = {a_debris - a_ref:.6e}  "
              f"({(a_debris / a_ref - 1) * 100:+.1f}% over the reference)")
        print(f"resistance = (a_x - a_x_ref) * mean(totalFluidMass) = {resistance:.6e}  (lattice force units)")

    else:  # fixed
        q_debris, q_debris_std, n_debris = steady_mean(args.debris, args.flow_probe, args.steady_from)
        q_ref, q_ref_std, n_ref = steady_mean(args.reference, args.flow_probe, args.steady_from)

        ratio = q_debris / q_ref if q_ref != 0 else float("nan")

        print(f"debris run    : {args.flow_probe} = {q_debris:.6e} +- {q_debris_std:.2e}  (n={n_debris} samples from t>={args.steady_from})")
        print(f"reference run : {args.flow_probe} = {q_ref:.6e} +- {q_ref_std:.2e}  (n={n_ref} samples from t>={args.steady_from})")
        print()
        print(f"Q/Q_ref = {ratio:.4f}  ({(ratio - 1) * 100:+.1f}% relative to the no-debris reference)")


if __name__ == "__main__":
    main()
