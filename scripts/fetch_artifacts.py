"""Download every artifact of one workflow run into results/<job>/.

    set GH_TOKEN=...            (PowerShell: $env:GH_TOKEN = "github_pat_...")
    python scripts/fetch_artifacts.py --run 44456107387

Each artifact holds m.out / p.out for one job id.
"""
import argparse
import io
import json
import os
import zipfile

import requests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
API = "https://api.github.com/repos/3264088435-code/auxetic-transferability"
HEAD = {"Accept": "application/vnd.github+json",
        "X-GitHub-Api-Version": "2022-11-28",
        "User-Agent": "doubao-qe"}


def get(url, token):
    # requests 会在跨主机重定向时丢掉 Authorization，artifact 又必须走一次
    # 签名 URL 的 302，所以这里不能用裸 urllib（会把头带过去，回 401）
    r = requests.get(url, headers={**HEAD, "Authorization": f"Bearer {token}"},
                     timeout=300, allow_redirects=True)
    r.raise_for_status()
    return r.content


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    args = ap.parse_args()
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if not token:
        raise SystemExit("GH_TOKEN not set")

    out_root = os.path.join(ROOT, "results")
    os.makedirs(out_root, exist_ok=True)

    arts, page = [], 1
    while True:
        data = json.loads(get(f"{API}/actions/runs/{args.run}/artifacts?per_page=100&page={page}", token))
        arts += data["artifacts"]
        if len(data["artifacts"]) < 100:
            break
        page += 1

    print(f"artifacts: {len(arts)}")
    for a in arts:
        name = a["name"]
        d = os.path.join(out_root, name)
        if os.path.exists(os.path.join(d, "p.out")) and os.path.exists(os.path.join(d, "m.out")):
            continue
        blob = get(f"{API}/actions/artifacts/{a['id']}/zip", token)
        os.makedirs(d, exist_ok=True)
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            z.extractall(d)
        print(f"  {name}: {os.listdir(d)}")


if __name__ == "__main__":
    main()
