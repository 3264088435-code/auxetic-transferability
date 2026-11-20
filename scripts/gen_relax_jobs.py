"""生成 InP relaxed-ion C44 的两个单侧作业（layer D，k=12，剪切 eps_yz，h=1%）。

与 A/B/C 层的唯一差别：calculation='relax'（固定晶胞、离子弛豫），
因此 sigma_xy 含内应变（Born-Kleinman）贡献，C44 是 relaxed-ion 值。

用法：python scripts/gen_relax_jobs.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_inputs as G

JOB_BASE = "D_inp_k12_relax_j3"
K = 12
DEG = 0.020
H = 0.01
MODE = 3          # eps_yz -> C44
CIF = "structures/mp-20351.cif"


def main():
    cell, sites = G.load_cif(CIF)
    species = sorted({s for s, _ in sites}, key=lambda x: x)
    newjobs = []
    for s, tag in ((-1, "m"), (1, "p")):
        jid = "%s_%s" % (JOB_BASE, tag)
        d = os.path.join(G.ROOT, "runs", jid)
        os.makedirs(d, exist_ok=True)
        G.write_input(os.path.join(d, tag + ".in"), prefix=jid, dft="PBE", k=K,
                      degauss=DEG, cell=G.deformed(cell, G.strain(MODE, s, H)),
                      sites=sites, species=species, calc="relax")
        with open(os.path.join(d, "pseudo.list"), "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(G.UPF[x] for x in species) + "\n")
        newjobs.append(dict(job=jid, layer="D", material="InP", k=K, degauss=DEG,
                            dft="PBE", mode=MODE, h=H, cell="prim", files=[tag + ".in"]))
        print("生成", jid, "nat=%d ntyp=%d" % (len(sites), len(species)))

    p = os.path.join(G.ROOT, "plan", "jobs.json")
    jobs = [j for j in json.load(open(p, encoding="utf-8")) if j["layer"] != "D"]
    jobs += newjobs
    jobs.sort(key=lambda x: x["job"])
    with open(p, "w", encoding="utf-8") as f:
        json.dump(jobs, f, indent=1, ensure_ascii=False)
    print("jobs.json 共", len(jobs), "条")


if __name__ == "__main__":
    main()
