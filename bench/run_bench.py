#!/usr/bin/env python3
"""Run remersher over the benchmark cases and score it against reference outputs.

For each case in cases.py this runs the remersher CLI, computes quality metrics for its output
and, when bench/reference/<mesh>__<paramset>.obj exists, for the reference output too. Results
land in bench/results/<timestamp>/ (results.json, results.csv, report.md) and are mirrored to
bench/results/latest/.

usage: run_bench.py [--bin build/remersher] [--jobs N] [--only bumpy__P0 ...]
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import datetime as dt
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import cases  # noqa: E402
import metrics  # noqa: E402

# Differences smaller than these (in the metric's own unit) count as a tie.
TIE_EPS = {
    "target_err_pct": 1.0, "nonquad_pct": 0.5, "irregular_pct": 0.5, "angle_dev_mean": 0.5,
    "angle_dev_p95": 1.0, "dev_mean_pct": 0.02, "hausdorff_pct": 0.1, "sharp_dev_pct": 0.02,
    "folded_pct": 0.1, "nonmanifold_edges": 0,
}
TIE_REL = 0.05

TABLE_COLS = [
    ("faces", "faces", "{:.0f}"), ("target_err_pct", "count err %", "{:.1f}"),
    ("nonquad_pct", "non-quad %", "{:.1f}"), ("irregular_pct", "poles %", "{:.1f}"),
    ("angle_dev_mean", "angle dev°", "{:.1f}"), ("dev_mean_pct", "mean dev %", "{:.3f}"),
    ("hausdorff_pct", "hausdorff %", "{:.2f}"), ("sharp_dev_pct", "crease dev %", "{:.3f}"),
    ("folded_pct", "folded %", "{:.2f}"),
]


def run_case(binary: str, mesh: str, ps: str, out_dir: str) -> dict:
    params = cases.PARAMSETS[ps]
    in_path = os.path.join(ROOT, "meshes", mesh + ".obj")
    inp = metrics.load_obj(in_path)
    T = inp.triangles()
    import numpy as np
    area = float(0.5 * np.linalg.norm(
        np.cross(inp.V[T[:, 1]] - inp.V[T[:, 0]], inp.V[T[:, 2]] - inp.V[T[:, 0]]), axis=1).sum())
    target = cases.target_count(params, len(T), area)
    args, unsupported = cases.to_cli_args(params, target)

    case = f"{mesh}__{ps}"
    out_path = os.path.join(out_dir, "meshes", case + ".obj")
    row = {"case": case, "mesh": mesh, "paramset": ps, "target": target,
           "unsupported": "; ".join(unsupported)}

    t0 = time.monotonic()
    try:
        proc = subprocess.run([binary, "-i", in_path, "-o", out_path, "-q", *args],
                              capture_output=True, text=True, timeout=600)
        row["seconds"] = round(time.monotonic() - t0, 2)
        row["error"] = proc.stderr.strip()[-300:] if proc.returncode else ""
    except subprocess.TimeoutExpired:
        row["seconds"] = 600.0
        row["error"] = "timeout"

    row["ours"] = None
    if not row["error"]:
        row["ours"] = metrics.evaluate(metrics.load_obj(out_path), inp, target)

    ref_path = os.path.join(HERE, "reference", case + ".obj")
    row["ref"] = metrics.evaluate(metrics.load_obj(ref_path), inp, target) if os.path.exists(ref_path) else None
    return row


def compare(ours: float, ref: float, key: str) -> int:
    """+1 if ours is better, -1 if worse, 0 for a tie (lower is better for all scored metrics)."""
    if any(isinstance(v, float) and math.isnan(v) for v in (ours, ref)):
        return 0
    diff = ref - ours
    if abs(diff) <= TIE_EPS[key] or abs(diff) <= TIE_REL * max(abs(ours), abs(ref)):
        return 0
    return 1 if diff > 0 else -1


def git_rev() -> str:
    try:
        rev = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                             capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=ROOT,
                               capture_output=True, text=True).stdout.strip()
        return (rev or "uncommitted") + ("+dirty" if dirty else "")
    except OSError:
        return "unknown"


def fmt(fmt_str: str, v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "–"
    return fmt_str.format(v)


def write_report(rows: list[dict], out_dir: str, rev: str) -> str:
    scored = [r for r in rows if r["ours"] and r["ref"]]
    lines = [f"# Remersher benchmark", "",
             f"- date: {dt.datetime.now():%Y-%m-%d %H:%M}", f"- remersher: `{rev}`",
             f"- cases: {len(rows)} run, {sum(1 for r in rows if r['error'])} failed, "
             f"{len(scored)} compared against a reference", ""]

    if scored:
        lines += ["## Head to head (lower is better)", "",
                  "| metric | wins | ties | losses | median ours | median ref |",
                  "|---|---:|---:|---:|---:|---:|"]
        total = [0, 0, 0]
        for key in metrics.LOWER_IS_BETTER:
            pairs = [(r["ours"].get(key), r["ref"].get(key)) for r in scored]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None
                     and not (isinstance(a, float) and math.isnan(a))
                     and not (isinstance(b, float) and math.isnan(b))]
            if not pairs:
                continue
            res = [compare(a, b, key) for a, b in pairs]
            w, t, l = res.count(1), res.count(0), res.count(-1)
            total = [total[0] + w, total[1] + t, total[2] + l]
            lines.append(f"| {key} | {w} | {t} | {l} | {statistics.median(a for a, _ in pairs):.3f} "
                         f"| {statistics.median(b for _, b in pairs):.3f} |")
        lines += [f"| **all metrics** | **{total[0]}** | **{total[1]}** | **{total[2]}** | | |", ""]

    lines += ["## Per case", "", "Each cell is `ours / reference`.", "",
              "| case | target | time s | " + " | ".join(c[1] for c in TABLE_COLS) + " | notes |",
              "|---|---:|---:|" + "---:|" * len(TABLE_COLS) + "---|"]
    for r in rows:
        cells = []
        for key, _, f in TABLE_COLS:
            a = r["ours"].get(key) if r["ours"] else None
            b = r["ref"].get(key) if r["ref"] else None
            mark = ""
            if a is not None and b is not None and key in TIE_EPS:
                mark = {1: " ✅", -1: " ❌", 0: ""}[compare(a, b, key)]
            cells.append(f"{fmt(f, a)} / {fmt(f, b)}{mark}")
        notes = r["error"] or (("unsupported: " + r["unsupported"]) if r["unsupported"] else "")
        lines.append(f"| {r['case']} | {r['target']} | {r['seconds']} | " + " | ".join(cells)
                     + f" | {notes} |")
    report = "\n".join(lines) + "\n"
    with open(os.path.join(out_dir, "report.md"), "w") as fh:
        fh.write(report)
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bin", default=os.path.join(ROOT, "build", "remersher"))
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2))
    ap.add_argument("--only", nargs="*", help="case names such as bumpy__P0")
    args = ap.parse_args()
    if not os.path.isfile(args.bin):
        sys.exit(f"remersher binary not found at {args.bin}; build it first")

    todo = cases.matrix()
    if args.only:
        todo = [c for c in todo if f"{c[0]}__{c[1]}" in args.only]

    stamp = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    out_dir = os.path.join(HERE, "results", stamp)
    os.makedirs(os.path.join(out_dir, "meshes"), exist_ok=True)

    rows = []
    with cf.ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(run_case, args.bin, m, p, out_dir): (m, p) for m, p in todo}
        for fut in cf.as_completed(futures):
            m, p = futures[fut]
            try:
                row = fut.result()
            except Exception as exc:  # metric failure should not sink the whole run
                row = {"case": f"{m}__{p}", "mesh": m, "paramset": p, "target": None,
                       "seconds": None, "unsupported": "", "error": f"bench error: {exc}",
                       "ours": None, "ref": None}
            rows.append(row)
            status = row["error"] or f"{row['ours']['faces']} faces"
            print(f"[{len(rows)}/{len(todo)}] {row['case']}: {status}", flush=True)

    order = {c: i for i, c in enumerate(f"{m}__{p}" for m, p in todo)}
    rows.sort(key=lambda r: order[r["case"]])
    rev = git_rev()

    with open(os.path.join(out_dir, "results.json"), "w") as fh:
        json.dump({"rev": rev, "rows": rows}, fh, indent=1)
    keys = sorted({k for r in rows for side in ("ours", "ref") if r[side] for k in r[side]})
    with open(os.path.join(out_dir, "results.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "target", "seconds", "error"] + [f"ours_{k}" for k in keys] + [f"ref_{k}" for k in keys])
        for r in rows:
            w.writerow([r["case"], r["target"], r["seconds"], r["error"]]
                       + [(r["ours"] or {}).get(k, "") for k in keys]
                       + [(r["ref"] or {}).get(k, "") for k in keys])

    report = write_report(rows, out_dir, rev)
    latest = os.path.join(HERE, "results", "latest")
    shutil.rmtree(latest, ignore_errors=True)
    shutil.copytree(out_dir, latest, ignore=shutil.ignore_patterns("meshes"))
    print("\n" + report)
    print(f"results: {out_dir}")


if __name__ == "__main__":
    main()
