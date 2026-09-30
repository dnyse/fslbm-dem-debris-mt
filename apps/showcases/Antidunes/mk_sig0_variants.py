#!/usr/bin/env python3
"""Make surface-tension-OFF (PhysicsParameters.surfaceTensionOverride 0) variants of the pinned-box and Bogner-like box scans (THESIS_PLAN.md,
2026-09-28: surface tension identified as the cause of the waterline excess buoyancy). Writes prms/*_sig0_gas.prm and jobs/run_72_*_sig0.sbatch;
the only prm difference is the one surfaceTensionOverride line."""
import re, glob
def sig0_prm(src):
    t = open(src, encoding="utf-8").read()
    assert "surfaceTensionOverride" not in t
    t, n = re.subn(r"(gravityOverride[^\n]*\n)", lambda m: m.group(1) + "   surfaceTensionOverride 0.0;   // surface-tension OFF test (nominal lattice value 0.000849592)\n", t, count=1); assert n == 1
    dst = src.replace("_gas.prm", "_sig0_gas.prm"); t = t.replace(src.split("/")[-1], dst.split("/")[-1], 1)
    open(dst, "w", encoding="utf-8").write(t); return dst
def sig0_sbatch(src, name_old, name_new, prm_re, lab_re, out_old, out_new, bk_old, bk_new):
    t = open(src, encoding="utf-8").read()
    t = re.sub(prm_re, r"\1_sig0_gas.prm", t); t = re.sub(lab_re, r"\1_sig0", t)
    for o, n in ((name_old, name_new), (out_old, out_new), (bk_old, bk_new)): assert o in t, o; t = t.replace(o, n)
    dst = src.replace(".sbatch", "_sig0.sbatch"); open(dst, "w", encoding="utf-8").write(t); return dst
n1 = [sig0_prm(f) for f in sorted(glob.glob("prms/pinbox_r*_gas.prm")) if "_sig0" not in f]
d1 = sig0_sbatch("jobs/run_72_pinbox_scan.sbatch", "antidunes-pinbox-scan", "antidunes-pinbox-scan-sig0", r"(pinbox_r\d+_a[\w.]+?)_gas\.prm", r"(box_r\d+_a[\w]+?)(?=\s*$)", "pinbox_scan_", "pinbox_scan_sig0_", '_pinbox"', '_pinbox_sig0"')
n2 = [sig0_prm(f) for f in sorted(glob.glob("prms/bognerlike_r*_gas.prm")) if "_sig0" not in f]
d2 = sig0_sbatch("jobs/run_72_bognerlike_box.sbatch", "antidunes-bognerlike-box", "antidunes-bognerlike-box-sig0", r"(bognerlike_r\d+_h\d+)_gas\.prm", r"(bl_r\d+_h\d+)(?=\s*$)", "bognerlike_box_", "bognerlike_box_sig0_", '_bognerlike"', '_bognerlike_sig0"')
print(len(n1), "pinbox prms;", d1); print(len(n2), "bognerlike prms;", d2)
