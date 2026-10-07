#!/usr/bin/env python3
"""Copy reference outputs from another remesher into bench/reference/.

Expects one folder per case named `<mesh>__<paramset>` containing an OBJ (default name
`out.geo.obj`), e.g. the output of running a licensed QuadRemesher over meshes/ with the
settings in cases.py. Reference outputs stay local: bench/reference/ is git-ignored.

usage: import_reference.py <runs_dir> [--name out.geo.obj]
"""
import argparse
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs_dir")
    ap.add_argument("--name", default="out.geo.obj", help="OBJ file name inside each case folder")
    args = ap.parse_args()

    dest = os.path.join(HERE, "reference")
    os.makedirs(dest, exist_ok=True)
    copied = 0
    for case in sorted(os.listdir(args.runs_dir)):
        src = os.path.join(args.runs_dir, case, args.name)
        if "__" in case and os.path.isfile(src):
            shutil.copyfile(src, os.path.join(dest, case + ".obj"))
            copied += 1
    print(f"copied {copied} reference meshes into {dest}")


if __name__ == "__main__":
    main()
