"""Quality metrics for quad remeshing results.

Every metric is computed only from the output mesh and the input mesh, so any remesher's
output can be scored the same way. Distances are reported as a percentage of the input's
bounding-box diagonal so results are comparable across meshes.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree


@dataclass
class PolyMesh:
    V: np.ndarray            # (n, 3) float64
    faces: list[list[int]]   # 0-based polygons

    def triangles(self) -> np.ndarray:
        tris = []
        for f in self.faces:
            tris.extend(triangulate_polygon(self.V, f))
        return np.asarray(tris, dtype=np.int64).reshape(-1, 3)


def triangulate_polygon(V: np.ndarray, poly: list[int]) -> list[tuple[int, int, int]]:
    """Ear clipping in the polygon's best-fit plane (same method as the engine's cleanup), so
    concave ngons are measured correctly; falls back to a fan for degenerate polygons."""
    n = len(poly)
    if n == 3:
        return [tuple(poly)]
    if n == 4:
        # Split along the diagonal that keeps both triangles facing the same way.
        P = V[poly]
        nrm = np.cross(P[2] - P[0], P[3] - P[1])
        t1 = np.cross(P[1] - P[0], P[2] - P[0]) @ nrm
        t2 = np.cross(P[2] - P[0], P[3] - P[0]) @ nrm
        if t1 > 0 and t2 > 0:
            return [(poly[0], poly[1], poly[2]), (poly[0], poly[2], poly[3])]
        return [(poly[1], poly[2], poly[3]), (poly[1], poly[3], poly[0])]
    P = V[poly]
    nrm = np.zeros(3)
    for i in range(n):
        a, b = P[i], P[(i + 1) % n]
        nrm += [(a[1] - b[1]) * (a[2] + b[2]), (a[2] - b[2]) * (a[0] + b[0]), (a[0] - b[0]) * (a[1] + b[1])]
    drop = int(np.argmax(np.abs(nrm)))
    ax, ay = (drop + 1) % 3, (drop + 2) % 3
    sign = 1.0 if nrm[drop] >= 0 else -1.0
    p2 = [(P[i][ax], sign * P[i][ay]) for i in range(n)]

    def cross2(a, b, c):
        return (p2[b][0] - p2[a][0]) * (p2[c][1] - p2[a][1]) - (p2[b][1] - p2[a][1]) * (p2[c][0] - p2[a][0])

    ring = list(range(n))
    tris = []
    while len(ring) > 3:
        m = len(ring)
        for k in range(m):
            a, b, c = ring[k - 1], ring[k], ring[(k + 1) % m]
            if cross2(a, b, c) <= 0:
                continue
            if any(q not in (a, b, c) and cross2(a, b, q) >= 0 and cross2(b, c, q) >= 0
                   and cross2(c, a, q) >= 0 for q in ring):
                continue
            tris.append((poly[a], poly[b], poly[c]))
            ring.pop(k)
            break
        else:
            return [(poly[0], poly[k], poly[k + 1]) for k in range(1, n - 1)]
    tris.append((poly[ring[0]], poly[ring[1]], poly[ring[2]]))
    return tris


def load_obj(path: str) -> PolyMesh:
    verts, faces = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                parts = line.split()
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith("f "):
                face = []
                for tok in line.split()[1:]:
                    i = int(tok.split("/")[0])
                    face.append(i - 1 if i > 0 else len(verts) + i)
                if len(face) >= 3:
                    faces.append(face)
    if not verts or not faces:
        raise ValueError(f"{path}: empty mesh")
    return weld(PolyMesh(np.asarray(verts, dtype=np.float64), faces))


def weld(mesh: PolyMesh, rel_tol: float = 1e-7) -> PolyMesh:
    """Merge vertices at (nearly) the same position. Some exporters (e.g. FBX round trips)
    split every face apart, which would otherwise make all topology metrics meaningless."""
    V = mesh.V
    diag = float(np.linalg.norm(V.max(0) - V.min(0))) or 1.0
    keys = np.round(V / (diag * rel_tol)).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    inverse = inverse.ravel()
    faces = []
    for f in mesh.faces:
        g = []
        for i in f:
            j = int(inverse[i])
            if not g or g[-1] != j:
                g.append(j)
        if len(g) > 1 and g[0] == g[-1]:
            g.pop()
        if len(g) >= 3:
            faces.append(g)
    return PolyMesh(V[first], faces)


# ---------------------------------------------------------------------------------------------
# geometry helpers

def _tri_normals(V, T):
    n = np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]])
    area2 = np.linalg.norm(n, axis=1)
    return n / np.maximum(area2, 1e-300)[:, None], 0.5 * area2


def sample_surface(V, T, count, rng):
    _, areas = _tri_normals(V, T)
    if areas.sum() <= 0:
        return V[T[:, 0]][:1].copy(), np.zeros(1, dtype=np.int64)
    idx = rng.choice(len(T), size=count, p=areas / areas.sum())
    r1 = np.sqrt(rng.random(count))[:, None]
    r2 = rng.random(count)[:, None]
    a, b, c = V[T[idx, 0]], V[T[idx, 1]], V[T[idx, 2]]
    return (1 - r1) * a + r1 * (1 - r2) * b + r1 * r2 * c, idx


def _closest_on_triangles(P, A, B, C):
    """Vectorised closest point on triangles ABC to points P (Ericson, RTCD 5.1.5)."""
    ab, ac, ap = B - A, C - A, P - A
    d1, d2 = (ab * ap).sum(-1), (ac * ap).sum(-1)
    bp = P - B
    d3, d4 = (ab * bp).sum(-1), (ac * bp).sum(-1)
    cp = P - C
    d5, d6 = (ab * cp).sum(-1), (ac * cp).sum(-1)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2

    denom = va + vb + vc
    denom = np.where(np.abs(denom) < 1e-300, 1e-300, denom)
    v = vb / denom
    w = vc / denom
    res = A + ab * v[..., None] + ac * w[..., None]

    def put(mask, value):
        nonlocal res
        res = np.where(mask[..., None], value, res)

    # Ericson tests regions in the order A, B, AB, C, AC, BC, interior and returns on the first
    # match; applying them in reverse lets earlier tests win, matching that priority.
    with np.errstate(divide="ignore", invalid="ignore"):
        m = (va <= 0) & ((d4 - d3) >= 0) & ((d5 - d6) >= 0)
        put(m, B + (C - B) * np.nan_to_num((d4 - d3) / ((d4 - d3) + (d5 - d6)))[..., None])
        m = (vb <= 0) & (d2 >= 0) & (d6 <= 0)
        put(m, A + ac * np.nan_to_num(d2 / (d2 - d6))[..., None])
        put((d6 >= 0) & (d5 <= d6), C)
        m = (vc <= 0) & (d1 >= 0) & (d3 <= 0)
        put(m, A + ab * np.nan_to_num(d1 / (d1 - d3))[..., None])
    put((d3 >= 0) & (d4 <= d3), B)
    put((d1 <= 0) & (d2 <= 0), A)
    return res


class SurfaceQuery:
    """Approximate-exact point-to-surface distance: exact distance to the k triangles whose
    centroids are nearest, after densifying large triangles via extra seed points."""

    def __init__(self, V, T, k=16):
        self.V, self.T, self.k = V, T, min(k, len(T))
        self.normals, _ = _tri_normals(V, T)
        cent = V[T].mean(axis=1)
        # Seed points: centroids plus vertices, so long thin triangles are still found.
        seeds = np.concatenate([cent, V[T[:, 0]], V[T[:, 1]], V[T[:, 2]]])
        self.seed_tri = np.concatenate([np.arange(len(T))] * 4)
        self.tree = cKDTree(seeds)

    def query(self, P, chunk=20000):
        dist = np.empty(len(P))
        tri = np.empty(len(P), dtype=np.int64)
        for s in range(0, len(P), chunk):
            p = P[s:s + chunk]
            _, nn = self.tree.query(p, k=self.k)
            cand = self.seed_tri[nn]                         # (m, k)
            A, B, C = (self.V[self.T[cand, j]] for j in range(3))
            q = _closest_on_triangles(p[:, None, :], A, B, C)
            d = np.linalg.norm(q - p[:, None, :], axis=-1)
            best = d.argmin(axis=1)
            dist[s:s + chunk] = d[np.arange(len(p)), best]
            tri[s:s + chunk] = cand[np.arange(len(p)), best]
        return dist, tri


def _edges(faces):
    """Map undirected edge -> number of incident faces."""
    count: dict[tuple[int, int], int] = {}
    for f in faces:
        n = len(f)
        for i in range(n):
            a, b = f[i], f[(i + 1) % n]
            if a == b:
                continue
            e = (a, b) if a < b else (b, a)
            count[e] = count.get(e, 0) + 1
    return count


def sharp_edge_samples(mesh: PolyMesh, angle_deg: float, spacing: float):
    """Points sampled along input edges whose dihedral angle exceeds angle_deg."""
    T = mesh.triangles()
    normals, _ = _tri_normals(mesh.V, T)
    edge_faces: dict[tuple[int, int], list[int]] = {}
    for fi, t in enumerate(T):
        for i in range(3):
            a, b = int(t[i]), int(t[(i + 1) % 3])
            edge_faces.setdefault((min(a, b), max(a, b)), []).append(fi)
    cos_thr = math.cos(math.radians(angle_deg))
    pts = []
    for (a, b), fs in edge_faces.items():
        if len(fs) != 2 or float(normals[fs[0]] @ normals[fs[1]]) > cos_thr:
            continue
        pa, pb = mesh.V[a], mesh.V[b]
        n = max(1, int(np.linalg.norm(pb - pa) / spacing))
        t = (np.arange(n) + 0.5) / n
        pts.append(pa + (pb - pa) * t[:, None])
    return np.concatenate(pts) if pts else np.zeros((0, 3))


def edge_samples(mesh: PolyMesh, spacing: float):
    pts = []
    for a, b in _edges(mesh.faces):
        pa, pb = mesh.V[a], mesh.V[b]
        n = max(1, int(np.linalg.norm(pb - pa) / spacing)) + 1
        t = np.linspace(0, 1, n)
        pts.append(pa + (pb - pa) * t[:, None])
    return np.concatenate(pts)


# ---------------------------------------------------------------------------------------------

# Metrics where a lower value is better; used by the benchmark to decide wins and losses.
LOWER_IS_BETTER = [
    "target_err_pct", "nonquad_pct", "irregular_pct", "angle_dev_mean", "angle_dev_p95",
    "dev_mean_pct", "hausdorff_pct", "sharp_dev_pct", "flipped_pct", "nonmanifold_edges",
]


def evaluate(output: PolyMesh, reference_input: PolyMesh, target: int | None = None,
             sharp_angle: float = 35.0, samples: int = 60000, seed: int = 1) -> dict:
    rng = np.random.default_rng(seed)
    Vin = reference_input.V
    Tin = reference_input.triangles()
    diag = float(np.linalg.norm(Vin.max(0) - Vin.min(0)))
    m: dict = {}

    # --- topology --------------------------------------------------------------------------
    sizes = np.array([len(f) for f in output.faces])
    m["faces"] = int(len(sizes))
    m["quads"] = int((sizes == 4).sum())
    m["tris"] = int((sizes == 3).sum())
    m["ngons"] = int((sizes > 4).sum())
    m["nonquad_pct"] = 100.0 * (1 - m["quads"] / max(1, m["faces"]))
    if target:
        m["target_err_pct"] = 100.0 * abs(m["faces"] - target) / target

    edges = _edges(output.faces)
    m["nonmanifold_edges"] = int(sum(1 for c in edges.values() if c > 2))
    valence = np.zeros(len(output.V), dtype=np.int64)
    on_border = np.zeros(len(output.V), dtype=bool)
    for (a, b), c in edges.items():
        valence[a] += 1
        valence[b] += 1
        if c == 1:
            on_border[a] = on_border[b] = True
    used = valence > 0
    regular = np.where(on_border, valence == 3, valence == 4)
    # Border corners (valence 2) are expected on open meshes and not counted as poles.
    regular |= on_border & (valence == 2)
    m["irregular_pct"] = 100.0 * float((~regular & used).sum()) / max(1, int(used.sum()))

    # --- quad shape -------------------------------------------------------------------------
    quads = np.array([f for f in output.faces if len(f) == 4], dtype=np.int64).reshape(-1, 4)
    if len(quads):
        P = output.V[quads]
        prev, nxt = np.roll(P, 1, axis=1) - P, np.roll(P, -1, axis=1) - P
        cosang = (prev * nxt).sum(-1) / np.maximum(
            np.linalg.norm(prev, axis=-1) * np.linalg.norm(nxt, axis=-1), 1e-300)
        dev = np.abs(np.degrees(np.arccos(np.clip(cosang, -1, 1))) - 90.0).ravel()
        m["angle_dev_mean"] = float(dev.mean())
        m["angle_dev_p95"] = float(np.percentile(dev, 95))
        lens = np.linalg.norm(nxt, axis=-1).ravel()
        m["edge_len_cv"] = float(lens.std() / max(lens.mean(), 1e-300))
    else:
        m["angle_dev_mean"] = m["angle_dev_p95"] = m["edge_len_cv"] = float("nan")

    # --- surface fidelity -------------------------------------------------------------------
    Tout = output.triangles()
    in_query, out_query = SurfaceQuery(Vin, Tin), SurfaceQuery(output.V, Tout)
    p_out, tri_out = sample_surface(output.V, Tout, samples, rng)
    p_in, _ = sample_surface(Vin, Tin, samples, rng)
    d_out_to_in, nearest_in = in_query.query(p_out)
    d_in_to_out, _ = out_query.query(p_in)
    m["dev_mean_pct"] = 100.0 * float(np.concatenate([d_out_to_in, d_in_to_out]).mean()) / diag
    m["hausdorff_pct"] = 100.0 * float(max(d_out_to_in.max(), d_in_to_out.max())) / diag

    out_normals, _ = _tri_normals(output.V, Tout)
    dots = (out_normals[tri_out] * in_query.normals[nearest_in]).sum(-1)
    m["flipped_pct"] = 100.0 * float((dots < 0).mean())

    # --- hard edges -------------------------------------------------------------------------
    spacing = diag / 3000.0
    sharp = sharp_edge_samples(reference_input, sharp_angle, spacing)
    m["sharp_samples"] = int(len(sharp))
    if len(sharp):
        d, _ = cKDTree(edge_samples(output, spacing)).query(sharp)
        m["sharp_dev_pct"] = 100.0 * float(d.mean()) / diag
    else:
        m["sharp_dev_pct"] = float("nan")
    return m
