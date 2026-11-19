"""Build the GitHub Actions job matrix from plan/jobs.json.

Prints two lines for $GITHUB_OUTPUT:

    matrix={"job": ["A_hfrh_k20_d010_pbe_j0", ...]}
    count=6
"""
import argparse
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan", default="A")
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    jobs = json.load(open(os.path.join(ROOT, "plan", "jobs.json"), encoding="utf-8"))
    plan = args.plan.strip().upper()
    if plan != "ALL":
        jobs = [j for j in jobs if j["layer"].upper() == plan]
    if args.only:
        jobs = [j for j in jobs if args.only in j["job"]]

    names = sorted(j["job"] for j in jobs)
    print(f"matrix={json.dumps({'job': names}, separators=(',', ':'))}")
    print(f"count={len(names)}")


if __name__ == "__main__":
    main()
