"""Benchmark cases: which meshes run with which settings.

Settings are written in QuadRemesher's public UI vocabulary so the same case can be run in both
tools. `to_cli_args` translates them to remersher options; settings remersher does not support
yet are reported as `unsupported` so the report can flag those cases.
"""
from __future__ import annotations

PARAMSETS = {
    "P0": dict(TargetQuadCount=2000, CurvatureAdaptivness=50, AutoDetectHardEdges=1),
    "P1": dict(TargetQuadCount=500, CurvatureAdaptivness=50, AutoDetectHardEdges=1),
    "P2": dict(TargetQuadCount=8000, CurvatureAdaptivness=50, AutoDetectHardEdges=1),
    "P3": dict(TargetQuadCount=2000, CurvatureAdaptivness=50, AutoDetectHardEdges=1, ExactQuadCount=1),
    "P4": dict(TargetQuadCount=2000, CurvatureAdaptivness=10, AutoDetectHardEdges=1),
    "P5": dict(TargetQuadCount=2000, CurvatureAdaptivness=90, AutoDetectHardEdges=1),
    "P6": dict(TargetEdgeLength=0.25, CurvatureAdaptivness=50, AutoDetectHardEdges=1),
    "P7": dict(TargetQuadCount=2000, CurvatureAdaptivness=50, AutoDetectHardEdges=0),
    "P8": dict(TargetQuadCount=2000, CurvatureAdaptivness=50, AutoDetectHardEdges=1, SymAxis="X"),
    "P9": dict(TargetQuadCountAsInputPercentage=20, CurvatureAdaptivness=50, AutoDetectHardEdges=1),
}

CORE = ["P0", "P1", "P2"]
EXTRA = ["P3", "P4", "P5", "P6", "P7", "P8", "P9"]
ALL_MESHES = ["sphere_4160", "sphere_32768", "cube", "ngon_cube", "stair", "torus",
              "plane_hole", "gear", "bumpy", "twospheres"]
EXTRA_MESHES = ["sphere_4160", "gear", "bumpy", "plane_hole"]


def matrix():
    cases = [(m, p) for m in ALL_MESHES for p in CORE]
    cases += [(m, p) for m in EXTRA_MESHES for p in EXTRA]
    return cases


def target_count(params: dict, input_faces: int, input_area: float) -> int:
    if "TargetQuadCount" in params:
        return params["TargetQuadCount"]
    if "TargetEdgeLength" in params:
        return max(1, round(input_area / params["TargetEdgeLength"] ** 2))
    if "TargetQuadCountAsInputPercentage" in params:
        return max(1, round(input_faces * params["TargetQuadCountAsInputPercentage"] / 100))
    raise ValueError(f"no target in {params}")


def to_cli_args(params: dict, target: int) -> tuple[list[str], list[str]]:
    args = ["--target", str(target)]
    unsupported = []
    if params.get("ExactQuadCount"):
        args += ["--passes", "5", "--tolerance", "0.005"]
    # remersher's adaptivity is on/off for now; low adaptiveness maps to uniform sizing.
    if params.get("CurvatureAdaptivness", 50) < 25:
        args.append("--no-adaptive")
    elif params.get("CurvatureAdaptivness", 50) not in (50,):
        unsupported.append(f"CurvatureAdaptivness={params['CurvatureAdaptivness']} (binary only)")
    if not params.get("AutoDetectHardEdges", 1):
        args.append("--no-hard-edges")
    if "SymAxis" in params:
        args += ["--symmetry", params["SymAxis"].lower()]
    return args, unsupported
