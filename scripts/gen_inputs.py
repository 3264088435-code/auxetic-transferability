"""Generate QE SCF inputs for the A/B/C layers.

Convention (identical to 03_程序/qe_k24_strain.ipynb, so new numbers are directly
comparable with the archived ones):

    eps[j,j]   = s*h            for j = 0,1,2      (uniaxial, Voigt eps_j = s*h)
    eps[k,l]   = 0.5*s*h        for j = 3,4,5      (shear,   Voigt eps_j = s*h)
    cell' = cell @ (I + eps),  s = +1 / -1,  h = 0.01
    C[:, j] = (sigma(+) - sigma(-)) / (2h),  sigma = -0.1 * QE stress (kbar -> GPa)

One job = one strain mode = the two amplitudes (s = -1, +1).  Frozen ions
(calculation='scf'), exactly as the functional/relax notebooks did.

Usage:
    python scripts/gen_inputs.py
"""
import json
import os

from pymatgen.core import Element, Structure

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
H = 0.01
ECUTWFC = 60.0
ECUTRHO = 480.0
DEGAUSS = 0.020
CONV_THR = "1.0d-8"
MIXING_BETA = 0.4

UPF = {
    "Hf": "Hf.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Rh": "Rh.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "In": "In.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "P": "P.pbe-n-kjpaw_psl.1.0.0.UPF",
    "B": "B.pbe-n-kjpaw_psl.1.0.0.UPF",
    "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Sn": "Sn.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Al": "Al.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Zr": "Zr.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Cu": "Cu.pbe-dn-kjpaw_psl.1.0.0.UPF",
}

# material -> (cif, strain modes).  Cubic needs only eps_xx (C11, C12) and eps_yz (C44).
MATERIALS = {
    "InP":   ("structures/mp-20351.cif", [0, 3]),
    "BN":    ("structures/mp-984.cif", [0, 1, 2, 3, 4, 5]),
    "SnS":   ("structures/mp-559676.cif", [0, 1, 2, 3, 4, 5]),
    "HfAl":  ("structures/mp-249.cif", [0, 1, 2, 3, 4, 5]),
    "Zr2Cu": ("structures/mp-193.cif", [0, 1, 2, 3, 4, 5]),
}

# A layer: HfRh, the one already published in the archive.  Base cell and
# fractional coordinates are copied verbatim from qe_2026-11-17_k24/.
HFRH = dict(
    name="HfRh",
    cell=[[3.249047, 0.0, 0.0],
          [0.0, 3.249047, 0.0],
          [0.0, 0.0, 3.249047]],
    sites=[("Hf", (0.0, 0.0, 0.0)), ("Rh", (0.5, 0.5, 0.5))],
    modes=[0, 3],
)
A_SETTINGS = [
    dict(k=20, degauss=0.010, dft="PBE"),
    dict(k=24, degauss=0.010, dft="PBE"),
    dict(k=24, degauss=0.020, dft="PBEsol"),
    # 存档那批用的是 h=1e-3；C' 只剩十几 GPa 时 1% 步长的三阶项不能白送，
    # 这一组专门用来核对步长影响（其余设置与上一行完全相同）
    dict(k=24, degauss=0.010, dft="PBE", h=0.001),
]

# 几何控制：vc-relax 得到本管线自己的平衡胞，再在该胞上测 C'
RELAX_CELL = os.path.join(ROOT, "structures", "hfrh_relaxed_cell.json")


def strain(j, s, h=H):
    eps = [[0.0] * 3 for _ in range(3)]
    if j < 3:
        eps[j][j] = s * h
    else:
        k, l = [(1, 2), (0, 2), (0, 1)][j - 3]
        eps[k][l] = eps[l][k] = 0.5 * s * h
    return eps


def deformed(cell, eps):
    return [[sum(cell[i][m] * ((1.0 if m == n else 0.0) + eps[m][n]) for m in range(3))
             for n in range(3)] for i in range(3)]


def write_input(path, *, prefix, dft, k, degauss, cell, sites, species, calc="scf"):
    lines = ["&CONTROL",
             f" calculation='{calc}'"]
    if calc != "scf":
        lines += [" restart_mode='from_scratch'",
                  " etot_conv_thr=1.0d-5",
                  " forc_conv_thr=1.0d-4"]
    lines += [f" prefix='{prefix}'",
             " outdir='./tmp'",
             " pseudo_dir='./pseudo'",
             " tstress=.true.",
             " tprnfor=.true.",
             " disk_io='none'",
             " verbosity='low'",
             "/",
             "&SYSTEM",
             " ibrav=0",
             f" nat={len(sites)}",
             f" ntyp={len(species)}",
             f" ecutwfc={ECUTWFC:.1f}",
             f" ecutrho={ECUTRHO:.1f}",
             f" input_dft='{dft}'",
             f" occupations='smearing', smearing='mv', degauss={degauss:.3f}",
             "/",
             "&ELECTRONS",
             f" conv_thr={CONV_THR}",
             f" mixing_beta={MIXING_BETA}",
             "/"]
    if calc != "scf":
        lines += ["&IONS", " ion_dynamics='bfgs'", "/",
                  "&CELL", " cell_dynamics='bfgs'", " press_conv_thr=0.1", "/"]
    lines += ["ATOMIC_SPECIES"]
    for sym in species:
        lines.append(f" {sym:<3} {Element(sym).atomic_mass:.3f} {UPF[sym]}")
    lines.append("CELL_PARAMETERS angstrom")
    for row in cell:
        lines.append(" " + " ".join(f"{v:.9f}" for v in row))
    lines.append("ATOMIC_POSITIONS crystal")
    for sym, (x, y, z) in sites:
        lines.append(f" {sym:<3} {x:.9f} {y:.9f} {z:.9f}")
    lines.append("K_POINTS automatic")
    lines.append(f" {k} {k} {k} 0 0 0")
    lines.append("")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))


def emit(jobs, job_id, layer, meta, cell, sites, species, modes, out_root, h=H, cell_tag="prim"):
    d = os.path.join(out_root, "runs", job_id)
    os.makedirs(d, exist_ok=True)
    for j in modes:
        for s, tag in ((-1, "m"), (1, "p")):
            write_input(os.path.join(d, f"{tag}.in"),
                        prefix=job_id, dft=meta["dft"], k=meta["k"],
                        degauss=meta["degauss"], cell=deformed(cell, strain(j, s, h)),
                        sites=sites, species=species)
        jobs.append(dict(job=job_id, layer=layer, material=meta["name"], k=meta["k"],
                         degauss=meta["degauss"], dft=meta["dft"], mode=j, h=h, cell=cell_tag,
                         files=["m.in", "p.in"]))
    plist = [UPF[s] for s in species]
    with open(os.path.join(d, "pseudo.list"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(plist) + "\n")


def load_cif(rel):
    st = Structure.from_file(os.path.join(ROOT, rel))
    return st.lattice.matrix.tolist(), [
        (site.species.elements[0].symbol, tuple(round(c, 9) for c in site.frac_coords))
        for site in st
    ]


def main():
    jobs = []
    # ---- A layer ----
    a_species = ["Hf", "Rh"]
    for st in A_SETTINGS:
        tag = "pbe" if st["dft"] == "PBE" else st["dft"].lower()
        h = st.get("h", H)
        htag = "" if abs(h - H) < 1e-12 else f"_h{int(round(h * 1000)):03d}"
        for j in HFRH["modes"]:
            jid = f"A_hfrh_k{st['k']}{htag}_d{int(round(st['degauss'] * 1000)):03d}_{tag}_j{j}"
            meta = dict(name=HFRH["name"], k=st["k"], degauss=st["degauss"], dft=st["dft"])
            emit(jobs, jid, "A", meta, HFRH["cell"], HFRH["sites"], a_species,
                 [j], ROOT, h=h)

    # ---- A 层几何控制 ----
    wdir = os.path.join(ROOT, "runs", "A_hfrh_k24_relax")
    os.makedirs(wdir, exist_ok=True)
    write_input(os.path.join(wdir, "r.in"), prefix="A_hfrh_k24_relax", dft="PBE", k=24,
                degauss=0.020, cell=HFRH["cell"], sites=HFRH["sites"],
                species=a_species, calc="vc-relax")
    with open(os.path.join(wdir, "pseudo.list"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(UPF[s] for s in a_species) + "\n")
    jobs.append(dict(job="A_hfrh_k24_relax", layer="A", material="HfRh", k=24,
                     degauss=0.020, dft="PBE", mode=-1, h=0.0, files=["r.in"]))

    if os.path.exists(RELAX_CELL):
        cell = json.load(open(RELAX_CELL, encoding="utf-8"))["cell"]
        meta = dict(name="HfRh", k=24, degauss=0.020, dft="PBE")
        for j in HFRH["modes"]:
            jid = f"A_hfrh_k24_rlx_d020_pbe_j{j}"
            emit(jobs, jid, "A", meta, cell, HFRH["sites"], a_species, [j], ROOT, h=0.001,
                 cell_tag="relaxed")

    # ---- B and C layers ----
    for layer, grids in (("B", (6, 12, 16)), ("C", (24,))):
        for name, (cif, modes) in MATERIALS.items():
            cell, sites = load_cif(cif)
            species = sorted({s for s, _ in sites}, key=lambda x: x)
            for k in grids:
                for j in modes:
                    jid = f"{layer}_{name.lower()}_k{k}_d020_pbe_j{j}"
                    meta = dict(name=name, k=k, degauss=DEGAUSS, dft="PBE")
                    emit(jobs, jid, layer, meta, cell, sites, species, [j], ROOT)

    jobs.sort(key=lambda x: x["job"])
    with open(os.path.join(ROOT, "plan", "jobs.json"), "w", encoding="utf-8") as f:
        json.dump(jobs, f, indent=1, ensure_ascii=False)
    counts = {}
    for jb in jobs:
        counts[jb["layer"]] = counts.get(jb["layer"], 0) + 1
    print(f"job files: {len(jobs)}  (A={counts.get('A',0)} B={counts.get('B',0)} C={counts.get('C',0)})")
    print(f"SCF runs:  {2*len(jobs)}")


if __name__ == "__main__":
    main()
