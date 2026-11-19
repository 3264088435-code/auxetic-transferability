"""Turn the downloaded QE outputs into an elastic-constant table.

Layout expected: results/<job_id>/{m.out,p.out}

    python scripts/parse_results.py

Writes results/qe_summary.csv with one row per job:
    job, material, k, degauss, dft, mode, sigma_m[6], sigma_p[6], column C[:,j]
"""
import csv
import glob
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KBAR_TO_GPA = 0.1


def stress_block(path):
    """Return the last stress tensor (kbar) reported in a QE .out file."""
    text = open(path, encoding="utf-8", errors="replace").read()
    idx = text.rfind("total   stress")
    if idx < 0:
        return None
    tail = text[idx:].splitlines()
    rows = []
    for line in tail[1:]:
        m = re.match(r"\s*([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+"
                     r"([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)", line)
        if not m:
            if rows:
                break
            continue
        rows.append([float(v) for v in m.groups()])
        if len(rows) == 3:
            break
    if len(rows) != 3:
        return None
    m = rows
    sxx, syy, szz = m[0][0], m[1][1], m[2][2]
    syz = 0.5 * (m[1][2] + m[2][1])
    sxz = 0.5 * (m[0][2] + m[2][0])
    sxy = 0.5 * (m[0][1] + m[1][0])
    # QE stress is kbar with compression positive -> flip sign, go to GPa
    return [-KBAR_TO_GPA * v for v in (sxx, syy, szz, syz, sxz, sxy)]


def main():
    meta = {j["job"]: j for j in json.load(open(os.path.join(ROOT, "plan", "jobs.json"), encoding="utf-8"))}
    rows = []
    for d in sorted(glob.glob(os.path.join(ROOT, "results", "*"))):
        jid = os.path.basename(d)
        if jid not in meta:
            continue
        sm = stress_block(os.path.join(d, "m.out"))
        sp = stress_block(os.path.join(d, "p.out"))
        if sm is None or sp is None:
            print(f"!! {jid}: missing stress block")
            continue
        col = [(p - m) / (2 * 0.01) for m, p in zip(sm, sp)]
        rec = dict(meta[jid])
        rec.update({"sigma_m": " ".join(f"{v:.3f}" for v in sm),
                    "sigma_p": " ".join(f"{v:.3f}" for v in sp),
                    "C_col": " ".join(f"{v:.2f}" for v in col)})
        rows.append(rec)

    out = os.path.join(ROOT, "results", "qe_summary.csv")
    keys = ["job", "layer", "material", "k", "degauss", "dft", "mode",
            "sigma_m", "sigma_p", "C_col"]
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in keys})
    print(f"rows: {len(rows)} -> {out}")


if __name__ == "__main__":
    main()
