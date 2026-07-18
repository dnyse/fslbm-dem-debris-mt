#!/usr/bin/env python3
"""Summarize per-floater F_contact+lub magnitude from an Antidunes run's stdout log.

Turns thousands of [BoxForces]/[CylForces]/[FloaterForces] lines into one row per floater
UID: peak |F_contact|, the step it occurred at, the step range over which |F_contact| stays
above a noise floor (the "contact window"), whether it drops back below the floor afterward,
and whether any NaN/Inf ever showed up. Doesn't judge physics -- just replaces log-scrolling
with a table you can compare against the .prm header's documented expected outcome.

Usage:
    python3 summarize_contact_log.py <logfile> [--noise-floor 1e-4]
"""
import argparse
import math
import re
import sys
from collections import defaultdict

LINE_RE = re.compile(
    r"\[(?P<shape>BoxForces|CylForces|FloaterForces)\]"
    r"\s+uid=(?P<uid>\d+)\s+step=(?P<step>\d+)"
    r".*?F_contact\+lub=\(\s*(?P<x>[^,]+),\s*(?P<y>[^,]+),\s*(?P<z>[^)]+)\)"
)


def parse_component(raw):
    """Parse one F_contact component; returns (value_or_None, is_bad)."""
    try:
        v = float(raw.strip())
    except ValueError:
        return None, True
    if math.isnan(v) or math.isinf(v):
        return v, True
    return v, False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logfile")
    ap.add_argument("--noise-floor", type=float, default=1e-4,
                     help="magnitude below which F_contact is treated as 'no contact' (default: 1e-4)")
    args = ap.parse_args()

    # per uid: shape tag, list of (step, magnitude), bad-value flag
    per_uid = defaultdict(lambda: {"shape": None, "samples": [], "bad": False})

    with open(args.logfile, "r", errors="replace") as f:
        for line in f:
            m = LINE_RE.search(line)
            if not m:
                continue
            uid = int(m.group("uid"))
            step = int(m.group("step"))
            entry = per_uid[uid]
            entry["shape"] = m.group("shape")

            comps = []
            bad = False
            for key in ("x", "y", "z"):
                v, is_bad = parse_component(m.group(key))
                bad = bad or is_bad
                comps.append(v)

            if bad:
                entry["bad"] = True
                continue
            mag = math.sqrt(sum(c * c for c in comps))
            entry["samples"].append((step, mag))

    if not per_uid:
        print(f"No [BoxForces]/[CylForces]/[FloaterForces] lines found in {args.logfile}.")
        sys.exit(1)

    header = f"{'uid':>10}  {'shape':<13}  {'n':>6}  {'peak|F_contact|':>16}  {'peak_step':>10}  {'contact_window':>20}  {'returns_to_0':>13}  {'NaN/Inf':>8}"
    print(header)
    print("-" * len(header))

    any_contact = False
    for uid in sorted(per_uid):
        entry = per_uid[uid]
        samples = entry["samples"]
        if not samples:
            print(f"{uid:>10}  {entry['shape'] or '?':<13}  {0:>6}  {'--':>16}  {'--':>10}  {'--':>20}  {'--':>13}  {'YES' if entry['bad'] else 'no':>8}")
            continue

        peak_step, peak_mag = max(samples, key=lambda s: s[1])
        touching = [s for s in samples if s[1] > args.noise_floor]

        if touching:
            any_contact = True
            window = f"{touching[0][0]}-{touching[-1][0]}"
            last_step, _ = samples[-1]
            returns = "yes" if (last_step - touching[-1][0] > 0 and samples[-1][1] <= args.noise_floor) or touching[-1][0] != last_step else "n/a (ends in contact)"
        else:
            window = "none"
            returns = "n/a (never touched)"

        print(f"{uid:>10}  {entry['shape']:<13}  {len(samples):>6}  {peak_mag:>16.6g}  {peak_step:>10}  {window:>20}  {returns:>13}  {'YES' if entry['bad'] else 'no':>8}")

    print()
    if not any_contact:
        print("No floater ever exceeded the noise floor -- consistent with a gap-pass negative control.")
    else:
        print("At least one floater showed a contact window -- compare timing/duration against the .prm header's expected outcome.")


if __name__ == "__main__":
    main()
