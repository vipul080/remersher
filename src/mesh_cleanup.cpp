#include "mesh_cleanup.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>
#include <stdexcept>
#include <unordered_map>

namespace remersher {
namespace {

using Vec3 = std::array<double, 3>;

double dist2(const Vec3& a, const Vec3& b) {
    double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

struct CellHash {
    size_t operator()(const std::array<int64_t, 3>& c) const {
        return (size_t)(c[0] * 73856093) ^ (size_t)(c[1] * 19349663) ^ (size_t)(c[2] * 83492791);
    }
};

// Maps each input vertex to a welded vertex id. Uses a uniform grid with cell size = tolerance
// and checks the 27 neighbouring cells, so pairs straddling a cell border are still merged.
std::vector<int> weldVertices(const std::vector<Vec3>& in, double tol, std::vector<Vec3>& out) {
    std::vector<int> remap(in.size(), -1);
    if (tol <= 0) {
        out = in;
        for (size_t i = 0; i < in.size(); ++i) remap[i] = (int)i;
        return remap;
    }
    const double tol2 = tol * tol;
    std::unordered_map<std::array<int64_t, 3>, std::vector<int>, CellHash> grid;
    auto cellOf = [&](const Vec3& p) {
        return std::array<int64_t, 3>{(int64_t)std::floor(p[0] / tol), (int64_t)std::floor(p[1] / tol),
                                      (int64_t)std::floor(p[2] / tol)};
    };
    for (size_t i = 0; i < in.size(); ++i) {
        const auto c = cellOf(in[i]);
        int found = -1;
        for (int dx = -1; dx <= 1 && found < 0; ++dx)
            for (int dy = -1; dy <= 1 && found < 0; ++dy)
                for (int dz = -1; dz <= 1 && found < 0; ++dz) {
                    auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
                    if (it == grid.end()) continue;
                    for (int id : it->second)
                        if (dist2(out[id], in[i]) <= tol2) { found = id; break; }
                }
        if (found < 0) {
            found = (int)out.size();
            out.push_back(in[i]);
            grid[c].push_back(found);
        }
        remap[i] = found;
    }
    return remap;
}

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

// Splits a polygon (indices into pts) into triangles by ear clipping in the polygon's best-fit
// plane, so concave ngons (stars, L-shapes, CAD caps) are triangulated without triangles that
// stick out of the polygon or fold over. Falls back to a fan if the polygon is too degenerate
// for clipping (self-intersecting or collinear).
std::vector<std::array<int, 3>> triangulatePolygon(const std::vector<Vec3>& pts,
                                                   const std::vector<int>& poly) {
    const size_t n = poly.size();
    std::vector<std::array<int, 3>> tris;
    if (n == 3) return {{poly[0], poly[1], poly[2]}};

    // Newell normal and a 2D frame in its plane.
    Vec3 nrm{0, 0, 0};
    for (size_t i = 0; i < n; ++i) {
        const Vec3 &a = pts[poly[i]], &b = pts[poly[(i + 1) % n]];
        nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
        nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
        nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    const int drop = std::fabs(nrm[0]) > std::fabs(nrm[1])
                         ? (std::fabs(nrm[0]) > std::fabs(nrm[2]) ? 0 : 2)
                         : (std::fabs(nrm[1]) > std::fabs(nrm[2]) ? 1 : 2);
    const int ax = (drop + 1) % 3, ay = (drop + 2) % 3;
    const double sign = nrm[drop] >= 0 ? 1.0 : -1.0;  // keep the polygon counter-clockwise in 2D
    std::vector<std::array<double, 2>> p2(n);
    for (size_t i = 0; i < n; ++i) p2[i] = {pts[poly[i]][ax], sign * pts[poly[i]][ay]};

    auto cross2 = [&](size_t a, size_t b, size_t c) {
        return (p2[b][0] - p2[a][0]) * (p2[c][1] - p2[a][1]) - (p2[b][1] - p2[a][1]) * (p2[c][0] - p2[a][0]);
    };
    std::vector<size_t> ring(n);
    for (size_t i = 0; i < n; ++i) ring[i] = i;

    size_t guard = 0;
    while (ring.size() > 3 && guard++ < n * n) {
        bool clipped = false;
        const size_t m = ring.size();
        for (size_t k = 0; k < m; ++k) {
            const size_t a = ring[(k + m - 1) % m], b = ring[k], c = ring[(k + 1) % m];
            if (cross2(a, b, c) <= 0) continue;  // reflex or degenerate corner
            bool contains = false;
            for (size_t j = 0; j < m && !contains; ++j) {
                const size_t q = ring[j];
                if (q == a || q == b || q == c) continue;
                contains = cross2(a, b, q) >= 0 && cross2(b, c, q) >= 0 && cross2(c, a, q) >= 0;
            }
            if (contains) continue;
            tris.push_back({poly[a], poly[b], poly[c]});
            ring.erase(ring.begin() + k);
            clipped = true;
            break;
        }
        if (!clipped) break;
    }
    if (ring.size() == 3) {
        tris.push_back({poly[ring[0]], poly[ring[1]], poly[ring[2]]});
        return tris;
    }
    tris.clear();  // clipping got stuck: fall back to a fan
    for (size_t k = 1; k + 1 < n; ++k) tris.push_back({poly[0], poly[k], poly[k + 1]});
    return tris;
}

// True if triangle t traverses the undirected edge (a, b) as a -> b.
bool runsForward(const std::array<int, 3>& t, int a, int b) {
    for (int i = 0; i < 3; ++i)
        if (t[i] == a && t[(i + 1) % 3] == b) return true;
    return false;
}

}  // namespace

TriangleMesh cleanupForRemeshing(const Mesh& mesh, double weldTolerance, CleanupStats* stats) {
    CleanupStats st;
    TriangleMesh out;

    // --- weld ------------------------------------------------------------------------------
    Vec3 lo{INFINITY, INFINITY, INFINITY}, hi{-INFINITY, -INFINITY, -INFINITY};
    for (const auto& v : mesh.vertices)
        for (int j = 0; j < 3; ++j) {
            if (!std::isfinite(v[j])) throw std::runtime_error("input has non-finite vertex positions");
            lo[j] = std::min(lo[j], v[j]);
            hi[j] = std::max(hi[j], v[j]);
        }
    const double diag = std::sqrt(dist2(lo, hi));
    const std::vector<int> remap = weldVertices(mesh.vertices, weldTolerance * diag, out.vertices);
    st.weldedVertices = (int)(mesh.vertices.size() - out.vertices.size());

    // --- triangulate, drop degenerate and duplicate triangles -------------------------------
    std::unordered_map<uint64_t, std::vector<int>> seen;  // sorted (a,b) key -> triangles
    std::vector<int> poly;
    for (const auto& face : mesh.faces) {
        poly.clear();
        for (int v : face)
            if (poly.empty() || poly.back() != remap[v]) poly.push_back(remap[v]);
        while (poly.size() > 1 && poly.front() == poly.back()) poly.pop_back();
        if (poly.size() < 3) { ++st.droppedTriangles; continue; }
        for (const auto& t : triangulatePolygon(out.vertices, poly)) {
            if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) { ++st.droppedTriangles; continue; }
            std::array<int, 3> s = t;
            std::sort(s.begin(), s.end());
            auto& bucket = seen[edgeKey(s[0], s[1])];
            bool dup = false;
            for (int other : bucket) {
                std::array<int, 3> o = out.triangles[other];
                std::sort(o.begin(), o.end());
                if (o == s) { dup = true; break; }
            }
            if (dup) { ++st.droppedTriangles; continue; }
            bucket.push_back((int)out.triangles.size());
            out.triangles.push_back(t);
        }
    }
    if (out.triangles.empty()) throw std::runtime_error("input mesh has no non-degenerate faces");

    // --- edge -> incident triangles ---------------------------------------------------------
    const int nt = (int)out.triangles.size();
    std::unordered_map<uint64_t, std::vector<int>> edgeTris;
    edgeTris.reserve(nt * 2);
    for (int f = 0; f < nt; ++f)
        for (int i = 0; i < 3; ++i)
            edgeTris[edgeKey(out.triangles[f][i], out.triangles[f][(i + 1) % 3])].push_back(f);

    // --- consistent winding per connected piece (BFS across manifold edges) ------------------
    std::vector<int> component(nt, -1);
    std::vector<char> flipped(nt, 0);
    auto flip = [&](int f) {
        std::swap(out.triangles[f][1], out.triangles[f][2]);
        flipped[f] ^= 1;
    };
    std::vector<int> members;
    for (int seed = 0; seed < nt; ++seed) {
        if (component[seed] >= 0) continue;
        const int comp = st.components++;
        members.clear();
        bool closed = true;
        std::queue<int> queue;
        component[seed] = comp;
        queue.push(seed);
        while (!queue.empty()) {
            const int f = queue.front();
            queue.pop();
            members.push_back(f);
            for (int i = 0; i < 3; ++i) {
                const int a = out.triangles[f][i], b = out.triangles[f][(i + 1) % 3];
                const auto& tris = edgeTris[edgeKey(a, b)];
                if (tris.size() != 2) { closed = false; continue; }  // border or non-manifold
                const int g = tris[0] == f ? tris[1] : tris[0];
                if (component[g] >= 0) continue;  // conflicts on visited faces: non-orientable, keep
                // A consistent neighbour runs the shared edge as b -> a.
                if (runsForward(out.triangles[g], a, b)) flip(g);
                component[g] = comp;
                queue.push(g);
            }
        }

        bool flipAll;
        if (closed) {
            double volume = 0;  // signed volume * 6; positive when faces point outwards
            for (int f : members) {
                const Vec3 &a = out.vertices[out.triangles[f][0]], &b = out.vertices[out.triangles[f][1]],
                           &c = out.vertices[out.triangles[f][2]];
                volume += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                          a[2] * (b[0] * c[1] - b[1] * c[0]);
            }
            flipAll = volume < 0;
        } else {
            size_t changed = 0;
            for (int f : members) changed += flipped[f];
            flipAll = changed * 2 > members.size();
        }
        if (flipAll)
            for (int f : members) flip(f);
    }
    for (int f = 0; f < nt; ++f) st.flippedTriangles += flipped[f];

    // --- drop unreferenced vertices ---------------------------------------------------------
    std::vector<int> used(out.vertices.size(), -1);
    std::vector<Vec3> compact;
    for (auto& t : out.triangles)
        for (int& v : t) {
            if (used[v] < 0) {
                used[v] = (int)compact.size();
                compact.push_back(out.vertices[v]);
            }
            v = used[v];
        }
    out.vertices = std::move(compact);

    if (stats) *stats = st;
    return out;
}

int repairPolygonMesh(Mesh& mesh, double weldTolerance) {
    if (mesh.vertices.empty()) return 0;
    Vec3 lo{INFINITY, INFINITY, INFINITY}, hi{-INFINITY, -INFINITY, -INFINITY};
    for (const auto& v : mesh.vertices)
        for (int j = 0; j < 3; ++j) lo[j] = std::min(lo[j], v[j]), hi[j] = std::max(hi[j], v[j]);
    std::vector<Vec3> welded;
    const std::vector<int> remap = weldVertices(mesh.vertices, weldTolerance * std::sqrt(dist2(lo, hi)), welded);

    const size_t before = mesh.faces.size();
    std::vector<std::vector<int>> faces;
    std::unordered_map<uint64_t, std::vector<int>> bySmallestEdge;  // duplicate detection
    for (const auto& f : mesh.faces) {
        std::vector<int> g;
        for (int v : f)
            if (g.empty() || g.back() != remap[v]) g.push_back(remap[v]);
        while (g.size() > 1 && g.front() == g.back()) g.pop_back();
        if (g.size() < 3) continue;
        std::vector<int> sorted = g;
        std::sort(sorted.begin(), sorted.end());
        if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) continue;  // folded
        auto& bucket = bySmallestEdge[edgeKey(sorted[0], sorted[1])];
        bool dup = false;
        for (int other : bucket) {
            std::vector<int> o = faces[other];
            std::sort(o.begin(), o.end());
            if (o == sorted) { dup = true; break; }
        }
        if (dup) continue;
        bucket.push_back((int)faces.size());
        faces.push_back(std::move(g));
    }

    // Edges with more than two faces: drop the smallest faces until the edge is manifold.
    auto area = [&](const std::vector<int>& f) {
        Vec3 n{0, 0, 0};
        for (size_t i = 0; i < f.size(); ++i) {
            const Vec3 &a = welded[f[i]], &b = welded[f[(i + 1) % f.size()]];
            n[0] += a[1] * b[2] - a[2] * b[1], n[1] += a[2] * b[0] - a[0] * b[2], n[2] += a[0] * b[1] - a[1] * b[0];
        }
        return 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    };
    std::vector<char> keep(faces.size(), 1);
    for (int round = 0; round < 4; ++round) {
        std::unordered_map<uint64_t, std::vector<int>> edgeFaces;
        for (int fi = 0; fi < (int)faces.size(); ++fi) {
            if (!keep[fi]) continue;
            const auto& f = faces[fi];
            for (size_t i = 0; i < f.size(); ++i) edgeFaces[edgeKey(f[i], f[(i + 1) % f.size()])].push_back(fi);
        }
        bool changed = false;
        for (auto& [k, fs] : edgeFaces) {
            if (fs.size() <= 2) continue;
            std::sort(fs.begin(), fs.end(), [&](int a, int b) { return area(faces[a]) > area(faces[b]); });
            for (size_t i = 2; i < fs.size(); ++i) keep[fs[i]] = 0, changed = true;
        }
        if (!changed) break;
    }

    Mesh out;
    std::vector<int> used(welded.size(), -1);
    for (size_t fi = 0; fi < faces.size(); ++fi) {
        if (!keep[fi]) continue;
        for (int& v : faces[fi]) {
            if (used[v] < 0) {
                used[v] = (int)out.vertices.size();
                out.vertices.push_back(welded[v]);
            }
            v = used[v];
        }
        out.faces.push_back(std::move(faces[fi]));
    }
    const int removed = (int)(before - out.faces.size());
    mesh = std::move(out);
    return removed;
}

}  // namespace remersher
