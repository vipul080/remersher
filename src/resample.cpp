#include "resample.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace remersher {
namespace {

using Vec3 = std::array<double, 3>;
using Tri = std::array<int, 3>;

Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Vec3 operator*(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 normalized(const Vec3& a) {
    double n = norm(a);
    return n > 0 ? a * (1.0 / n) : Vec3{0, 0, 0};
}

uint64_t key(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

Vec3 triNormal(const Vec3& a, const Vec3& b, const Vec3& c) { return cross(b - a, c - a); }

// Closest point on triangle abc to p (Ericson, Real-Time Collision Detection 5.1.5).
Vec3 closestOnTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    Vec3 ab = b - a, ac = c - a, ap = p - a;
    double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return a;
    Vec3 bp = p - b;
    double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return b;
    double vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) return a + ab * (d1 / (d1 - d3));
    Vec3 cp = p - c;
    double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return c;
    double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) return a + ac * (d2 / (d2 - d6));
    double va = d3 * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));
    double denom = 1.0 / (va + vb + vc);
    return a + ab * (vb * denom) + ac * (vc * denom);
}

// Projects points onto the original surface. Triangles are binned into a uniform grid; a query
// looks at the point's cell and its 26 neighbours, which is enough because smoothing moves a
// vertex by less than one target edge length (the cell size is two of them).
class SurfaceProjector {
   public:
    SurfaceProjector(const std::vector<Vec3>& V, const std::vector<Tri>& F, double cell)
        : V_(V), F_(F), cell_(cell) {
        for (int f = 0; f < (int)F.size(); ++f) {
            Vec3 lo = V[F[f][0]], hi = lo;
            for (int k = 1; k < 3; ++k)
                for (int j = 0; j < 3; ++j)
                    lo[j] = std::min(lo[j], V[F[f][k]][j]), hi[j] = std::max(hi[j], V[F[f][k]][j]);
            auto c0 = cellOf(lo), c1 = cellOf(hi);
            double count = (double)(c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1);
            if (count > 4096) { large_.push_back(f); continue; }  // checked on every query
            for (int64_t x = c0[0]; x <= c1[0]; ++x)
                for (int64_t y = c0[1]; y <= c1[1]; ++y)
                    for (int64_t z = c0[2]; z <= c1[2]; ++z) cells_[hash({x, y, z})].push_back(f);
        }
    }

    bool project(const Vec3& p, Vec3& out) const {
        double best = INFINITY;
        auto consider = [&](int f) {
            Vec3 q = closestOnTriangle(p, V_[F_[f][0]], V_[F_[f][1]], V_[F_[f][2]]);
            Vec3 d = q - p;
            double d2 = dot(d, d);
            if (d2 < best) best = d2, out = q;
        };
        auto c = cellOf(p);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = cells_.find(hash({c[0] + dx, c[1] + dy, c[2] + dz}));
                    if (it != cells_.end())
                        for (int f : it->second) consider(f);
                }
        for (int f : large_) consider(f);
        return best < INFINITY;
    }

   private:
    std::array<int64_t, 3> cellOf(const Vec3& p) const {
        return {(int64_t)std::floor(p[0] / cell_), (int64_t)std::floor(p[1] / cell_),
                (int64_t)std::floor(p[2] / cell_)};
    }
    static uint64_t hash(const std::array<int64_t, 3>& c) {
        return ((uint64_t)(c[0] & 0x1FFFFF) << 42) | ((uint64_t)(c[1] & 0x1FFFFF) << 21) |
               (uint64_t)(c[2] & 0x1FFFFF);
    }
    const std::vector<Vec3>& V_;
    const std::vector<Tri>& F_;
    double cell_;
    std::unordered_map<uint64_t, std::vector<int>> cells_;
    std::vector<int> large_;
};

class Remesher {
   public:
    Remesher(const TriangleMesh& in, double target, double hardAngle)
        : V(in.vertices), F(in.triangles), alive(in.triangles.size(), 1), L(target),
          projector(in.vertices, in.triangles, 2 * target) {
        buildAdjacency();
        // Feature edges: borders, non-manifold edges and (optionally) creases.
        const double cosHard = hardAngle > 0 ? std::cos(hardAngle * M_PI_ / 180.0) : -2.0;
        for (const auto& [k, faces] : edgeFaces) {
            bool feature = faces.size() != 2;
            if (!feature && cosHard > -2.0) {
                Vec3 n0 = normalized(faceNormal(faces[0])), n1 = normalized(faceNormal(faces[1]));
                feature = dot(n0, n1) < cosHard;
            }
            if (feature) features.insert(k);
            if (faces.size() > 2) {
                frozen.insert((int)(k >> 32));
                frozen.insert((int)(k & 0xffffffff));
            }
        }
    }

    TriangleMesh run(int iterations, ResampleStats& st) {
        const double hi = 4.0 / 3.0 * L, lo = 4.0 / 5.0 * L;
        for (int it = 0; it < iterations; ++it) {
            for (int round = 0; round < 40; ++round) {
                int n = splitPass(hi);
                st.splits += n;
                if (!n) break;
            }
            for (int round = 0; round < 10; ++round) {
                int n = collapsePass(lo, hi);
                st.collapses += n;
                if (!n) break;
            }
            st.flips += flipPass();
            smoothPass();
            ++st.iterations;
        }
        return compact();
    }

   private:
    std::vector<Vec3> V;
    std::vector<Tri> F;
    std::vector<char> alive;
    double L;
    SurfaceProjector projector;

    std::unordered_map<uint64_t, std::vector<int>> edgeFaces;
    std::vector<std::vector<int>> vertFaces;
    std::unordered_set<uint64_t> features;
    std::unordered_set<int> frozen;  // vertices on non-manifold edges
    std::vector<char> touched;

    static constexpr double M_PI_ = 3.14159265358979323846;

    Vec3 faceNormal(int f) const { return triNormal(V[F[f][0]], V[F[f][1]], V[F[f][2]]); }

    void buildAdjacency() {
        edgeFaces.clear();
        vertFaces.assign(V.size(), {});
        for (int f = 0; f < (int)F.size(); ++f) {
            if (!alive[f]) continue;
            for (int i = 0; i < 3; ++i) {
                edgeFaces[key(F[f][i], F[f][(i + 1) % 3])].push_back(f);
                vertFaces[F[f][i]].push_back(f);
            }
        }
        touched.assign(V.size(), 0);
    }

    int featureDegree(int v) const {
        int n = 0;
        for (int u : neighbours(v)) n += features.count(key(v, u)) ? 1 : 0;
        return n;
    }

    bool locked(int v) const { return frozen.count(v) || featureDegree(v) > 0; }

    std::vector<int> neighbours(int v) const {
        std::vector<int> n;
        for (int f : vertFaces[v])
            for (int u : F[f])
                if (u != v) n.push_back(u);
        std::sort(n.begin(), n.end());
        n.erase(std::unique(n.begin(), n.end()), n.end());
        return n;
    }

    void touch(int v) {
        touched[v] = 1;
        for (int u : neighbours(v)) touched[u] = 1;
    }

    int splitPass(double hi) {
        buildAdjacency();
        std::vector<std::pair<double, uint64_t>> longEdges;
        for (const auto& [k, faces] : edgeFaces) {
            if (faces.size() > 2) continue;
            int a = (int)(k >> 32), b = (int)(k & 0xffffffff);
            double len = norm(V[a] - V[b]);
            if (len > hi) longEdges.push_back({len, k});
        }
        std::sort(longEdges.rbegin(), longEdges.rend());
        int count = 0;
        for (const auto& [len, k] : longEdges) {
            int a = (int)(k >> 32), b = (int)(k & 0xffffffff);
            const auto& faces = edgeFaces[k];
            bool busy = touched[a] || touched[b];
            for (int f : faces)
                for (int v : F[f]) busy |= touched[v] != 0;
            if (busy) continue;

            const int m = (int)V.size();
            V.push_back((V[a] + V[b]) * 0.5);
            vertFaces.emplace_back();
            touched.push_back(1);
            for (int f : faces) {
                // Rotate so the split edge is (x, y) in the face's winding order.
                Tri t = F[f];
                for (int r = 0; r < 3 && !((t[0] == a && t[1] == b) || (t[0] == b && t[1] == a)); ++r)
                    t = {t[1], t[2], t[0]};
                F[f] = {t[0], m, t[2]};
                F.push_back({m, t[1], t[2]});
                alive.push_back(1);
                touched[t[2]] = 1;
            }
            if (features.count(k)) {
                features.erase(k);
                features.insert(key(a, m));
                features.insert(key(m, b));
            }
            touched[a] = touched[b] = 1;
            ++count;
        }
        return count;
    }

    // Removes r by merging it into k. Returns false (changing nothing) if that would break the
    // topology, flip a face, create an over-long edge, or move a feature off its feature line.
    bool tryCollapse(int r, int k, uint64_t ek, double hi) {
        if (frozen.count(r)) return false;
        const int rFeat = featureDegree(r);
        if (rFeat > 0 && !(features.count(ek) && rFeat == 2)) return false;

        const auto& ef = edgeFaces[ek];
        if (ef.empty() || ef.size() > 2) return false;
        std::vector<int> opposite;
        for (int f : ef)
            for (int v : F[f])
                if (v != r && v != k) opposite.push_back(v);
        std::sort(opposite.begin(), opposite.end());

        // Link condition: r and k may only share the vertices opposite the collapsed edge.
        auto nr = neighbours(r), nk = neighbours(k);
        std::vector<int> common;
        std::set_intersection(nr.begin(), nr.end(), nk.begin(), nk.end(), std::back_inserter(common));
        if (common != opposite) return false;
        // Do not pinch a border: an interior edge between two border vertices.
        if (ef.size() == 2 && featureDegree(k) > 0 && rFeat > 0 && !features.count(ek)) return false;
        if (nr.size() + nk.size() > 20) return false;  // keep valences reasonable

        for (int n : nr)
            if (n != k && norm(V[n] - V[k]) > hi) return false;
        for (int f : vertFaces[r]) {
            if (F[f][0] == k || F[f][1] == k || F[f][2] == k) continue;
            Vec3 before = faceNormal(f);
            Tri t = F[f];
            for (int& v : t)
                if (v == r) v = k;
            Vec3 after = triNormal(V[t[0]], V[t[1]], V[t[2]]);
            double nb = norm(before), na = norm(after);
            if (na <= 1e-14 * L * L || dot(before, after) < 0.5 * nb * na) return false;
        }

        for (int n : nr) touched[n] = 1;
        for (int n : nk) touched[n] = 1;
        touched[r] = touched[k] = 1;
        for (int f : vertFaces[r]) {
            if (F[f][0] == k || F[f][1] == k || F[f][2] == k) {
                alive[f] = 0;
                continue;
            }
            for (int& v : F[f])
                if (v == r) v = k;
        }
        // Feature edges (r, x) become (k, x).
        for (int n : nr)
            if (n != k && features.count(key(r, n))) features.insert(key(k, n));
        return true;
    }

    int collapsePass(double lo, double hi) {
        buildAdjacency();
        std::vector<std::pair<double, uint64_t>> shortEdges;
        for (const auto& [k, faces] : edgeFaces) {
            int a = (int)(k >> 32), b = (int)(k & 0xffffffff);
            double len = norm(V[a] - V[b]);
            if (len < lo) shortEdges.push_back({len, k});
        }
        std::sort(shortEdges.begin(), shortEdges.end());
        int count = 0;
        for (const auto& [len, k] : shortEdges) {
            int a = (int)(k >> 32), b = (int)(k & 0xffffffff);
            if (touched[a] || touched[b]) continue;
            // Prefer removing the vertex that is free to move.
            bool aLocked = locked(a), bLocked = locked(b);
            if (aLocked && !bLocked) std::swap(a, b);
            if (tryCollapse(b, a, k, hi) || tryCollapse(a, b, k, hi)) ++count;
        }
        return count;
    }

    int flipPass() {
        buildAdjacency();
        std::vector<int> valence(V.size(), 0);
        std::vector<char> border(V.size(), 0);
        for (const auto& [k, faces] : edgeFaces) {
            int a = (int)(k >> 32), b = (int)(k & 0xffffffff);
            ++valence[a], ++valence[b];
            if (faces.size() != 2) border[a] = border[b] = 1;
        }
        auto target = [&](int v) { return border[v] ? 4 : 6; };
        int count = 0;
        std::vector<uint64_t> keys;
        keys.reserve(edgeFaces.size());
        for (const auto& kv : edgeFaces) keys.push_back(kv.first);
        std::sort(keys.begin(), keys.end());  // deterministic order
        for (uint64_t k : keys) {
            const auto& ef = edgeFaces[k];
            if (ef.size() != 2 || features.count(k)) continue;
            // Name the corners so face 0 runs a -> b -> c and face 1 runs b -> a -> d.
            const Tri& t0 = F[ef[0]];
            const Tri& t1 = F[ef[1]];
            int a = -1, b = -1, c = -1, d = -1;
            const int e0 = (int)(k >> 32), e1 = (int)(k & 0xffffffff);
            for (int i = 0; i < 3; ++i)
                if ((t0[i] == e0 && t0[(i + 1) % 3] == e1) || (t0[i] == e1 && t0[(i + 1) % 3] == e0))
                    a = t0[i], b = t0[(i + 1) % 3], c = t0[(i + 2) % 3];
            for (int j = 0; j < 3; ++j)
                if (t1[j] == b && t1[(j + 1) % 3] == a) d = t1[(j + 2) % 3];
            if (a < 0 || d < 0) continue;  // inconsistent winding: leave alone
            if (c == d || touched[a] || touched[b] || touched[c] || touched[d]) continue;
            if (edgeFaces.count(key(c, d))) continue;
            if (valence[a] <= 3 || valence[b] <= 3) continue;

            auto dev = [&](int v, int val) { return std::abs(val - target(v)); };
            int before = dev(a, valence[a]) + dev(b, valence[b]) + dev(c, valence[c]) + dev(d, valence[d]);
            int after = dev(a, valence[a] - 1) + dev(b, valence[b] - 1) + dev(c, valence[c] + 1) +
                        dev(d, valence[d] + 1);
            if (after >= before) continue;

            Vec3 n0 = faceNormal(ef[0]), n1 = faceNormal(ef[1]);
            if (dot(normalized(n0), normalized(n1)) < 0.9) continue;  // would change the shape
            Tri f0 = {d, b, c}, f1 = {c, a, d};
            Vec3 m0 = triNormal(V[f0[0]], V[f0[1]], V[f0[2]]), m1 = triNormal(V[f1[0]], V[f1[1]], V[f1[2]]);
            Vec3 avg = normalized(n0 + n1);
            if (dot(m0, avg) <= 1e-12 * L * L || dot(m1, avg) <= 1e-12 * L * L) continue;

            F[ef[0]] = f0;
            F[ef[1]] = f1;
            --valence[a], --valence[b], ++valence[c], ++valence[d];
            touched[a] = touched[b] = touched[c] = touched[d] = 1;
            ++count;
        }
        return count;
    }

    void smoothPass() {
        buildAdjacency();
        std::vector<Vec3> next = V;
        for (int v = 0; v < (int)V.size(); ++v) {
            if (vertFaces[v].empty() || locked(v)) continue;
            auto nb = neighbours(v);
            if (nb.empty()) continue;
            Vec3 centroid{0, 0, 0}, n{0, 0, 0};
            for (int u : nb) centroid = centroid + V[u];
            centroid = centroid * (1.0 / nb.size());
            for (int f : vertFaces[v]) n = n + faceNormal(f);
            n = normalized(n);
            Vec3 d = centroid - V[v];
            Vec3 p = V[v] + (d - n * dot(n, d));
            Vec3 q;
            next[v] = projector.project(p, q) ? q : p;
        }
        V = std::move(next);
    }

    TriangleMesh compact() const {
        TriangleMesh out;
        std::vector<int> remap(V.size(), -1);
        for (int f = 0; f < (int)F.size(); ++f) {
            if (!alive[f]) continue;
            Tri t = F[f];
            if (t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) continue;
            for (int& v : t) {
                if (remap[v] < 0) {
                    remap[v] = (int)out.vertices.size();
                    out.vertices.push_back(V[v]);
                }
                v = remap[v];
            }
            out.triangles.push_back(t);
        }
        return out;
    }
};

}  // namespace

double sliverFraction(const TriangleMesh& mesh) {
    if (mesh.triangles.empty()) return 0;
    const double cos15 = std::cos(15.0 * 3.14159265358979323846 / 180.0);
    size_t slivers = 0;
    for (const auto& t : mesh.triangles) {
        for (int i = 0; i < 3; ++i) {
            const Vec3& p = mesh.vertices[t[i]];
            Vec3 u = mesh.vertices[t[(i + 1) % 3]] - p, w = mesh.vertices[t[(i + 2) % 3]] - p;
            double nu = norm(u), nw = norm(w);
            if (nu == 0 || nw == 0 || dot(u, w) / (nu * nw) > cos15) {
                ++slivers;
                break;
            }
        }
    }
    return (double)slivers / mesh.triangles.size();
}

TriangleMesh resampleIsotropic(const TriangleMesh& mesh, double targetEdge, double hardEdgeAngle,
                               int iterations, ResampleStats* stats) {
    ResampleStats st;
    Remesher r(mesh, targetEdge, hardEdgeAngle);
    TriangleMesh out = r.run(iterations, st);
    if (stats) *stats = st;
    return out;
}

}  // namespace remersher
