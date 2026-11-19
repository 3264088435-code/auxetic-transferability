"""Assemble 6x6 elastic tensors from results/qe_summary.csv and judge them.

    python scripts/build_tensors.py

Conventions (same as 03_程序/qe_hfrh.ipynb):
  * column j of C comes from (sigma(+) - sigma(-)) / (2h), h = 0.01
  * C is symmetrised and projected onto the crystal's point group (pymatgen)
  * Born stability: every eigenvalue of C positive
  * nu_min: direction scan over a uniform grid in cos(theta) and phi

Cubic jobs only carry modes 0 (eps_xx) and 3 (eps_yz), so their tensor is rebuilt
from C11 / C12 / C44 before any projection.
"""
import csv
import json
import os

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.tensors import Tensor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIFS = {"InP": "structures/mp-20351.cif", "BN": "structures/mp-984.cif",
        "SnS": "structures/mp-559676.cif", "HfAl": "structures/mp-249.cif",
        "Zr2Cu": "structures/mp-193.cif"}
HFRH_CELL = np.array([[3.249047, 0, 0], [0, 3.249047, 0], [0, 0, 3.249047]])


def rotations_for(material):
    if material == "HfRh":
        from pymatgen.core import Lattice
        st = Structure(Lattice(HFRH_CELL),
                       ["Hf", "Rh"], [[0, 0, 0], [0.5, 0.5, 0.5]])
        sga = SpacegroupAnalyzer(st, symprec=1e-3, angle_tolerance=5)
    else:
        st = Structure.from_file(os.path.join(ROOT, CIFS[material]))
        sga = SpacegroupAnalyzer(st, symprec=1e-3, angle_tolerance=5)
    return [op.rotation_matrix for op in sga.get_point_group_operations(cartesian=True)]


def project(C6, rots):
    """Project the 6x6 tensor onto the crystal point group, pure numpy.

    C6 -> 4-index C_ijkl without factors (sigma_ij = C_ijkl eps_kl sums over all
    k,l, and eps_4 in Voigt is the doubled shear, which is why the Voigt matrix
    maps straight across).  Averaging over the point group removes the
    symmetry-breaking noise of the finite differences.
    """
    pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    Ct = np.zeros((3, 3, 3, 3))
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            v = C6[a, b]
            Ct[i, j, k, l] = Ct[j, i, k, l] = Ct[i, j, l, k] = Ct[j, i, l, k] = v
    acc = np.zeros_like(Ct)
    for R in rots:
        acc += np.einsum("ip,jq,kr,ls,pqrs->ijkl", R, R, R, R, Ct)
    Ct = acc / len(rots)
    return np.array([[Ct[i, j, k, l] for (k, l) in pairs] for (i, j) in pairs])


def cubic_from_columns(cols):
    C11 = cols[0][0]
    C12 = 0.5 * (cols[0][1] + cols[0][2])
    C44 = cols[3][3]
    C = np.zeros((6, 6))
    C[:3, :3] = C12
    np.fill_diagonal(C[:3, :3], C11)
    np.fill_diagonal(C[3:, 3:], C44)
    return C


def directions(n_theta=181, n_phi=361):
    ct = 1.0 - 2.0 * (np.arange(n_theta) + 0.5) / n_theta
    st = np.sqrt(np.clip(1 - ct * ct, 0, None))
    ph = 2 * np.pi * np.arange(n_phi) / n_phi
    n = np.stack([st[:, None] * np.cos(ph)[None, :],
                  st[:, None] * np.sin(ph)[None, :],
                  np.repeat(ct[:, None], n_phi, axis=1)], axis=-1).reshape(-1, 3)
    return n / np.linalg.norm(n, axis=1)[:, None]


def s_tensor(S6):
    """6x6 compliance -> 4-index S_ijkl.

    Unlike C, the compliance needs the half factors: the Voigt shear strain
    component is the doubled shear, so S_yzyz = S44/4 for a cubic crystal.
    """
    pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
    g = [1.0, 1.0, 1.0, 0.5, 0.5, 0.5]
    S = np.zeros((3, 3, 3, 3))
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            v = S6[a, b] * g[a] * g[b]
            S[i, j, k, l] = S[j, i, k, l] = S[i, j, l, k] = S[j, i, l, k] = v
    return S


def nu_extremes(C6):
    S = s_tensor(np.linalg.inv(C6))
    n = directions()
    sig = n[:, :, None] * n[:, None, :]
    eps = np.einsum("ijkl,nkl->nij", S, sig)
    e_nn = np.einsum("nij,ni,nj->n", eps, n, n)

    helper = np.where(np.abs(n[:, 2:3]) > 0.9, [[1.0, 0, 0]], [[0, 0, 1.0]])
    a = np.cross(n, helper)
    a /= np.linalg.norm(a, axis=1)[:, None]
    b = np.cross(n, a)
    m11 = np.einsum("ni,nij,nj->n", a, eps, a)
    m22 = np.einsum("ni,nij,nj->n", b, eps, b)
    m12 = np.einsum("ni,nij,nj->n", a, eps, b)
    half = 0.5 * (m11 + m22)
    disc = np.sqrt(np.clip(0.25 * (m11 - m22) ** 2 + m12 ** 2, 0, None))
    nu = np.concatenate([-(half + disc) / e_nn, -(half - disc) / e_nn])
    return float(nu.min()), float(nu.max())


def main():
    src = os.path.join(ROOT, "results", "qe_summary.csv")
    rows = list(csv.DictReader(open(src, encoding="utf-8")))

    groups = {}
    for r in rows:
        groups.setdefault((r["material"], int(r["k"]), r["dft"]), {})[int(r["mode"])] = \
            [float(v) for v in r["C_col"].split()]

    out = []
    for (mat, k, dft), cols in sorted(groups.items()):
        if mat in ("HfRh", "InP") and set(cols) == {0, 3}:
            C = cubic_from_columns(cols)
        else:
            C = np.zeros((6, 6))
            for j, v in cols.items():
                C[:, j] = v
                C[j, :] = v
        C = project(0.5 * (C + C.T), rotations_for(mat))
        C = 0.5 * (C + C.T)
        eig = np.linalg.eigvalsh(C)
        born = bool(eig.min() > 0)
        nu_min, nu_max = nu_extremes(C) if born else (float("nan"), float("nan"))
        rec = dict(material=mat, k=k, dft=dft, n_modes=len(cols),
                   min_eig=round(float(eig.min()), 2), born=born,
                   nu_min=round(nu_min, 4), nu_max=round(nu_max, 4),
                   C11=round(C[0, 0], 1), C22=round(C[1, 1], 1), C33=round(C[2, 2], 1),
                   C12=round(C[0, 1], 1), C13=round(C[0, 2], 1), C23=round(C[1, 2], 1),
                   C44=round(C[3, 3], 1), C55=round(C[4, 4], 1), C66=round(C[5, 5], 1))
        if mat == "HfRh":
            rec["Cp"] = round(0.5 * (C[0, 0] - C[0, 1]), 1)
            rec["cubic_thr"] = round(C[3, 3] * C[0, 0] / (C[0, 0] + 2 * C[0, 1]), 1)
        out.append(rec)
        print(json.dumps(rec, ensure_ascii=False))

    path = os.path.join(ROOT, "results", "tensors.csv")
    keys = ["material", "k", "dft", "n_modes", "min_eig", "born", "nu_min", "nu_max",
            "C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66"]
    extra = [k for k in ("Cp", "cubic_thr") if any(k in r for r in out)]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys + extra)
        w.writeheader()
        for r in out:
            w.writerow({k: r.get(k, "") for k in keys + extra})
    print(f"\n{len(out)} tensors -> {path}")


if __name__ == "__main__":
    main()
