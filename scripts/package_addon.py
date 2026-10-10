#!/usr/bin/env python3
"""Builds dist/remersher_blender-<platform>.zip: the Blender add-on with the remersher CLI bundled
in its bin/ folder. Install it in Blender via Preferences > Add-ons > Install from Disk.

usage: package_addon.py [path/to/remersher[.exe]]
"""
import os
import platform
import shutil
import stat
import sys
import tempfile
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def platform_tag():
    system = {"Darwin": "darwin", "Linux": "linux", "Windows": "windows"}.get(platform.system(), platform.system().lower())
    machine = platform.machine().lower()
    machine = {"amd64": "x86_64", "aarch64": "arm64"}.get(machine, machine)
    return f"{system}-{machine}"


def main():
    exe = "remersher.exe" if os.name == "nt" else "remersher"
    candidates = [sys.argv[1]] if len(sys.argv) > 1 else [
        os.path.join(ROOT, "build", exe), os.path.join(ROOT, "build", "Release", exe)]
    binary = next((c for c in candidates if os.path.isfile(c)), None)
    if not binary:
        sys.exit(f"remersher binary not found (looked in {candidates}); build it first")

    out_dir = os.path.join(ROOT, "dist")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"remersher_blender-{platform_tag()}.zip")
    with tempfile.TemporaryDirectory() as stage:
        pkg = os.path.join(stage, "remersher_blender")
        shutil.copytree(os.path.join(ROOT, "blender", "remersher_blender"), pkg,
                        ignore=shutil.ignore_patterns("__pycache__"))
        os.makedirs(os.path.join(pkg, "bin"), exist_ok=True)
        shutil.copy2(binary, os.path.join(pkg, "bin", exe))
        for name in ("LICENSE", "THIRD_PARTY_NOTICES.md"):
            shutil.copy2(os.path.join(ROOT, name), pkg)
        if os.path.exists(out):
            os.remove(out)
        with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
            for folder, _, files in os.walk(pkg):
                for f in files:
                    path = os.path.join(folder, f)
                    info = zipfile.ZipInfo.from_file(path, os.path.relpath(path, stage))
                    if f == exe:  # keep the executable bit for unzip tools that honour it
                        info.external_attr = (stat.S_IFREG | 0o755) << 16
                    with open(path, "rb") as fh:
                        zf.writestr(info, fh.read(), zipfile.ZIP_DEFLATED)
    print(out)


if __name__ == "__main__":
    main()
