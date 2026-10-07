#include "quality.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace remersher {
namespace {

using Vec3 = std::array<double, 3>;

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// Uniform hash grid for nearest-point queries.
class PointGrid {
   public:
    PointGrid(const std::vector<Vec3>& pts, double cell) : pts_(pts), cell_(cell) {
        for (int i = 0; i < (int)pts.size(); ++i) cells_[key(cellOf(pts[i]))].push_back(i);
    }

    // Index of the nearest point, or -1 if none lies within maxDist.
    int nearest(const Vec3& p, double maxDist, double* outDist = nullptr) const {
        const auto c = cellOf(p);
        const int rings = (int)std::ceil(maxDist / cell_);
        int best = -1;
        double best2 = maxDist * maxDist;
        for (int r = 0; r <= rings; ++r) {
            // Points in ring r are at least (r - 1) * cell away; stop once that exceeds the best.
            if (best >= 0 && (r - 1) * cell_ > std::sqrt(best2)) break;
            for (int dx = -r; dx <= r; ++dx)
                for (int dy = -r; dy <= r; ++dy)
                    for (int dz = -r; dz <= r; ++dz) {
                        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != r) continue;
                        auto it = cells_.find(key({c[0] + dx, c[1] + dy, c[2] + dz}));
                        if (it == cells_.end()) continue;
                        for (int i : it->second) {
                            Vec3 d = sub(pts_[i], p);
                            double d2 = dot(d, d);
                            if (d2 <= best2) { best2 = d2; best = i; }
                        }
                    }
        }
        if (outDist) *outDist = std::sqrt(best2);
        return best;
    }

   private:
    std::array<int64_t, 3> cellOf(const Vec3& p) const {
        return {(int64_t)std::floor(p[0] / cell_), (int64_t)std::floor(p[1] / cell_),
                (int64_t)std::floor(p[2] / cell_)};
    }
    static uint64_t key(const std::array<int64_t, 3>& c) {
        return ((uint64_t)(c[0] & 0x1FFFFF) << 42) | ((uint64_t)(c[1] & 0x1FFFFF) << 21) |
               (uint64_t)(c[2] & 0x1FFFFF);
    }

    const std::vector<Vec3>& pts_;
    double cell_;
    std::unordered_map<uint64_t, std::vector<int>> cells_;
};

}  // namespace

bool QualityCheck::acceptable(std::string* why) const {
    auto fail = [&](const std::string& msg) {
        if (why) *why = msg;
        return false;
    };
    if (farOutputFraction > 0.01) return fail("output strays from the surface");
    if (uncoveredFraction > 0.01) return fail("output leaves holes");
    if (flippedFraction > 0.03) return fail("output has folded faces");
    return true;
}

QualityCheck checkQuality(const TriangleMesh& input, const Mesh& output) {
    QualityCheck q;
    if (output.faces.empty() || input.triangles.empty()) {
        q.uncoveredFraction = 1;
        return q;
    }

    Vec3 lo{INFINITY, INFINITY, INFINITY}, hi{-INFINITY, -INFINITY, -INFINITY};
    for (const auto& v : input.vertices)
        for (int j = 0; j < 3; ++j) lo[j] = std::min(lo[j], v[j]), hi[j] = std::max(hi[j], v[j]);
    const double diag = norm(sub(hi, lo));

    // Mean output edge length (the expected spacing of output vertices).
    double outEdge = 0;
    size_t outEdges = 0;
    for (const auto& f : output.faces)
        for (size_t i = 0; i < f.size(); ++i, ++outEdges)
            outEdge += norm(sub(output.vertices[f[(i + 1) % f.size()]], output.vertices[f[i]]));
    outEdge /= std::max<size_t>(1, outEdges);

    // Sample the input surface uniformly by area; each sample keeps its triangle's normal. Input
    // vertices alone are not enough: coarse CAD-like inputs have large faces with no interior
    // vertices, where holes would go unnoticed and vertex normals blur across hard edges.
    std::vector<double> cumArea;
    std::vector<Vec3> faceNormal;
    cumArea.reserve(input.triangles.size());
    double total = 0;
    for (const auto& t : input.triangles) {
        const auto &a = input.vertices[t[0]], &b = input.vertices[t[1]], &c = input.vertices[t[2]];
        Vec3 n = cross(sub(b, a), sub(c, a));
        double len = norm(n);
        total += 0.5 * len;
        cumArea.push_back(total);
        faceNormal.push_back(len > 0 ? Vec3{n[0] / len, n[1] / len, n[2] / len} : Vec3{0, 0, 0});
    }
    const size_t nSamples = std::clamp<size_t>(8 * output.faces.size(), 20000, 200000);
    std::vector<Vec3> samples;
    std::vector<int> sampleTri;
    samples.reserve(nSamples);
    uint64_t state = 0x9E3779B97F4A7C15ull;  // fixed seed: the check must be deterministic
    auto rnd = [&]() {
        state ^= state << 13, state ^= state >> 7, state ^= state << 17;
        return (state >> 11) * (1.0 / 9007199254740992.0);
    };
    for (size_t i = 0; i < nSamples; ++i) {
        int t = (int)(std::lower_bound(cumArea.begin(), cumArea.end(), rnd() * total) - cumArea.begin());
        t = std::min(t, (int)input.triangles.size() - 1);
        double r1 = std::sqrt(rnd()), r2 = rnd();
        const auto& tri = input.triangles[t];
        const auto &a = input.vertices[tri[0]], &b = input.vertices[tri[1]], &c = input.vertices[tri[2]];
        Vec3 p;
        for (int j = 0; j < 3; ++j) p[j] = (1 - r1) * a[j] + r1 * (1 - r2) * b[j] + r1 * r2 * c[j];
        samples.push_back(p);
        sampleTri.push_back(t);
    }
    const double spacing = std::sqrt(total / nSamples);

    // Output vertices must sit on the input surface.
    const PointGrid inGrid(samples, std::max(spacing, diag * 1e-4));
    const double farLimit = 0.01 * diag + 2 * spacing;
    size_t far = 0;
    for (const auto& v : output.vertices)
        if (inGrid.nearest(v, farLimit) < 0) ++far;
    q.farOutputFraction = (double)far / output.vertices.size();

    // Every surface sample should be near some output vertex; a point inside a quad is at most
    // ~0.71 * outEdge from a corner, and 2.5x leaves room for adaptive sizing.
    const PointGrid outGrid(output.vertices, std::max(outEdge, diag * 1e-4));
    const double holeLimit = 2.5 * outEdge + spacing;
    size_t uncovered = 0;
    for (const auto& p : samples)
        if (outGrid.nearest(p, holeLimit) < 0) ++uncovered;
    q.uncoveredFraction = (double)uncovered / samples.size();

    // Face orientation against the input face under it, plus quad corner angles. Faces whose
    // normal is nearly perpendicular to the input (straddling a hard edge) are not counted.
    size_t flipped = 0, corners = 0;
    double angleSum = 0;
    for (const auto& f : output.faces) {
        Vec3 n{0, 0, 0}, center{0, 0, 0};
        for (size_t i = 0; i < f.size(); ++i) {
            const auto& a = output.vertices[f[i]];
            const auto& b = output.vertices[f[(i + 1) % f.size()]];
            Vec3 c = cross(a, b);
            for (int j = 0; j < 3; ++j) n[j] += c[j], center[j] += a[j] / f.size();
        }
        const double nlen = norm(n);
        const int nearest = inGrid.nearest(center, diag);
        if (nearest >= 0 && nlen > 0 && dot(n, faceNormal[sampleTri[nearest]]) < -0.3 * nlen) ++flipped;
        if (f.size() == 4) {
            for (size_t i = 0; i < 4; ++i, ++corners) {
                const auto& p = output.vertices[f[i]];
                Vec3 u = sub(output.vertices[f[(i + 3) % 4]], p), w = sub(output.vertices[f[(i + 1) % 4]], p);
                double denom = norm(u) * norm(w);
                double ang = denom > 0 ? std::acos(std::clamp(dot(u, w) / denom, -1.0, 1.0)) : 0;
                angleSum += std::abs(ang * 180.0 / 3.14159265358979323846 - 90.0);
            }
        }
    }
    q.flippedFraction = (double)flipped / output.faces.size();
    q.angleDeviation = corners ? angleSum / corners : 0;
    return q;
}

}  // namespace remersher
