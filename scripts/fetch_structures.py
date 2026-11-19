"""Fetch conventional-symmetry cells for the B/C layer materials from the Materials Project.

Run once, locally.  The CIFs it writes are committed to the repo so CI never needs an API key.

    python scripts/fetch_structures.py

Requires: pymatgen, an MP API key in D:\\aiwork\\cail\\.mp_key
"""
import json
import os
import urllib.parse
import urllib.request

from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KEYFILE = r"D:\aiwork\cail\.mp_key"

# old (pre-renumbering) MP ids that the 5,000-pair table uses
TARGETS = {
    "InP": ("mp-20351", "F-43m", 216),
    "BN": ("mp-984", "P6_3/mmc", 194),
    "SnS": ("mp-559676", "Cmcm", 63),
    "HfAl": ("mp-249", "Cmcm", 63),
    "Zr2Cu": ("mp-193", "I4/mmm", 139),
}


def fetch(ids):
    key = open(KEYFILE).read().strip()
    url = (
        "https://api.materialsproject.org/materials/summary/?"
        + urllib.parse.urlencode(
            {"material_ids": ",".join(ids),
             "_fields": "material_id,formula_pretty,symmetry,nsites,structure",
             "_limit": 50}
        )
    )
    req = urllib.request.Request(
        url,
        headers={
            "X-API-KEY": key,
            "accept": "application/json",
            "User-Agent": "Mozilla/5.0 (compatible; doubao-qe/1.0)",
        },
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())["data"]


def main():
    out_dir = os.path.join(ROOT, "structures")
    os.makedirs(out_dir, exist_ok=True)

    data = fetch([v[0] for v in TARGETS.values()])
    by_formula = {}
    for m in data:
        by_formula.setdefault(m["formula_pretty"], []).append(m)

    index = []
    for name, (mid, want, want_no) in TARGETS.items():
        cands = by_formula.get(name, [])
        # the pair table pins the space group; pick the entry whose symmetry matches
        hit = [c for c in cands if c["symmetry"]["symbol"] == want]
        if len(hit) != 1:
            raise SystemExit(f"{name}: expected exactly 1 {want}, got {len(hit)} of {len(cands)}")
        m = hit[0]

        struct = Structure.from_dict(m["structure"])
        conv = SpacegroupAnalyzer(struct, symprec=1e-3).get_conventional_standard_structure()
        conv = SpacegroupAnalyzer(conv, symprec=1e-3).get_conventional_standard_structure()

        cif = os.path.join(out_dir, f"{mid}.cif")
        CifWriter(conv, symprec=1e-3).write_file(cif)

        sga = SpacegroupAnalyzer(conv, symprec=1e-3)
        got_no = sga.get_space_group_number()
        if got_no != want_no:
            raise SystemExit(f"{name}: conventional cell fell back to #{got_no} "
                             f"{sga.get_space_group_symbol()}, expected #{want_no} {want}")
        lat = conv.lattice
        rec = {
            "name": name,
            "mp_id_table": mid,
            "mp_id_current": m["material_id"],
            "spacegroup": sga.get_space_group_symbol(),
            "number": sga.get_space_group_number(),
            "nsites": len(conv),
            "abc": [round(lat.a, 4), round(lat.b, 4), round(lat.c, 4)],
            "angles": [round(lat.alpha, 3), round(lat.beta, 3), round(lat.gamma, 3)],
            "cif": os.path.relpath(cif, ROOT).replace("\\", "/"),
        }
        index.append(rec)
        print(json.dumps(rec, ensure_ascii=False))

    with open(os.path.join(out_dir, "index.json"), "w", encoding="utf-8") as f:
        json.dump(index, f, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
