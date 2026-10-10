#include "crease.h"

#include <cmath>
#include <cstdint>
#include <unordered_map>

namespace remersher {
namespace {

using Vec3 = std::array<double, 3>;

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

}  // namespace

std::vector<std::array<Vec3, 2>> detectCreases(const TriangleMesh& mesh, double angleDegrees) {
    std::vector<std::array<Vec3, 2>> out;
    if (angleDegrees <= 0 || angleDegrees >= 180) return out;
    const auto& V = mesh.vertices;
    const auto& T = mesh.triangles;

    std::vector<Vec3> normal(T.size());
    for (size_t f = 0; f < T.size(); ++f) {
        const Vec3 &a = V[T[f][0]], &b = V[T[f][1]], &c = V[T[f][2]];
        const Vec3 u{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        Vec3 n{u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]};
        const double l = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        normal[f] = l > 0 ? Vec3{n[0] / l, n[1] / l, n[2] / l} : Vec3{0, 0, 0};
    }
    std::unordered_map<uint64_t, std::vector<int>> edgeFaces;
    for (int f = 0; f < (int)T.size(); ++f)
        for (int i = 0; i < 3; ++i) edgeFaces[edgeKey(T[f][i], T[f][(i + 1) % 3])].push_back(f);

    // Dihedral angle (degrees) of every manifold edge; borders and non-manifold edges get none.
    std::unordered_map<uint64_t, double> dihedral;
    dihedral.reserve(edgeFaces.size());
    for (const auto& [k, fs] : edgeFaces) {
        if (fs.size() != 2) continue;
        const Vec3 &n0 = normal[fs[0]], &n1 = normal[fs[1]];
        const double c = std::max(-1.0, std::min(1.0, n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2]));
        dihedral[k] = std::acos(c) * 180.0 / 3.14159265358979323846;
    }

    // Candidate edges by angle, then judged per connected group: real creases form lines and
    // closed loops (cube edges, the rim of a dent) with few loose ends, while a coarse smooth
    // surface yields scattered short pieces with many loose ends.
    std::vector<uint64_t> candidates;
    for (const auto& [k, theta] : dihedral)
        if (theta > angleDegrees) candidates.push_back(k);
    std::unordered_map<int, std::vector<int>> incident;  // vertex -> candidate indices
    for (int i = 0; i < (int)candidates.size(); ++i) {
        incident[(int)(candidates[i] >> 32)].push_back(i);
        incident[(int)(candidates[i] & 0xffffffff)].push_back(i);
    }
    std::vector<int> group(candidates.size(), -1);
    for (int seed = 0; seed < (int)candidates.size(); ++seed) {
        if (group[seed] >= 0) continue;
        std::vector<int> members{seed}, stack{seed};
        group[seed] = seed;
        while (!stack.empty()) {
            const int e = stack.back();
            stack.pop_back();
            for (int v : {(int)(candidates[e] >> 32), (int)(candidates[e] & 0xffffffff)})
                for (int other : incident[v])
                    if (group[other] < 0) group[other] = seed, stack.push_back(other), members.push_back(other);
        }
        std::unordered_map<int, int> degree;
        for (int e : members) ++degree[(int)(candidates[e] >> 32)], ++degree[(int)(candidates[e] & 0xffffffff)];
        int loose = 0;
        for (const auto& [v, d] : degree) loose += d == 1;
        if (members.size() < 4 || loose > 0.25 * degree.size()) continue;
        for (int e : members) out.push_back({V[candidates[e] >> 32], V[candidates[e] & 0xffffffff]});
    }
    return out;
}

}  // namespace remersher
