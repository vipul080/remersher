#!/usr/bin/env bash
# Builds dist/remersher_blender-<platform>.zip: the Blender add-on with the remersher CLI bundled
# in its bin/ folder. Install it in Blender via Preferences > Add-ons > Install from Disk.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
binary="${1:-$root/build/remersher}"
[ -x "$binary" ] || { echo "remersher binary not found at $binary; build it first" >&2; exit 1; }

platform="$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m)"
stage="$(mktemp -d)"
trap 'rm -rf "$stage"' EXIT
cp -R "$root/blender/remersher_blender" "$stage/"
find "$stage" -name __pycache__ -prune -exec rm -rf {} +
mkdir -p "$stage/remersher_blender/bin"
cp "$binary" "$stage/remersher_blender/bin/"
cp "$root/LICENSE" "$root/THIRD_PARTY_NOTICES.md" "$stage/remersher_blender/"

mkdir -p "$root/dist"
out="$root/dist/remersher_blender-$platform.zip"
rm -f "$out"
(cd "$stage" && zip -qr "$out" remersher_blender)
echo "$out"
