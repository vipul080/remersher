#include "symmetry.h"

#include <cmath>
#include <cstdint>
#include <unordered_map>

namespace remersher {
namespace {

uint64_t edgeKey(int a, int b) {
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint32_t)b;
}

}  // namespace

TriangleMesh clipToHalfSpace(const TriangleMesh& mesh, int axis) {
    TriangleMesh out;
    std::vector<int> keep(mesh.vertices.size(), -1);
    auto inside = [&](int v) { return mesh.vertices[v][axis] >= 0; };
    auto keptVertex = [&](int v) {
        if (keep[v] < 0) {
            keep[v] = (int)out.vertices.size();
            out.vertices.push_back(mesh.vertices[v]);
        }
        return keep[v];
    };
    std::unordered_map<uint64_t, int> cut;  // crossing edge -> vertex on the plane
    auto cutVertex = [&](int a, int b) {
        auto it = cut.find(edgeKey(a, b));
        if (it != cut.end()) return it->second;
        const auto &pa = mesh.vertices[a], &pb = mesh.vertices[b];
        const double t = pa[axis] / (pa[axis] - pb[axis]);
        std::array<double, 3> p;
        for (int j = 0; j < 3; ++j) p[j] = pa[j] + t * (pb[j] - pa[j]);
        p[axis] = 0.0;  // exactly on the plane, so the mirrored halves weld
        const int id = (int)out.vertices.size();
        out.vertices.push_back(p);
        cut[edgeKey(a, b)] = id;
        return id;
    };

    for (const auto& t : mesh.triangles) {
        // Sutherland-Hodgman against the half-space.
        std::vector<int> poly;
        for (int i = 0; i < 3; ++i) {
            const int a = t[i], b = t[(i + 1) % 3];
            const bool ia = inside(a), ib = inside(b);
            if (ia) poly.push_back(keptVertex(a));
            if (ia != ib && mesh.vertices[a][axis] != 0 && mesh.vertices[b][axis] != 0) poly.push_back(cutVertex(a, b));
        }
        for (size_t k = 1; k + 1 < poly.size(); ++k) {
            std::array<int, 3> tri{poly[0], poly[k], poly[k + 1]};
            if (tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]) out.triangles.push_back(tri);
        }
    }
    return out;
}

Mesh mirrorAcross(const Mesh& half, int axis, double seamTolerance) {
    Mesh out = half;
    const size_t n = half.vertices.size();

    // Border vertices near the plane become seam vertices, exactly on it.
    std::unordered_map<uint64_t, int> edgeCount;
    for (const auto& f : half.faces)
        for (size_t i = 0; i < f.size(); ++i) ++edgeCount[edgeKey(f[i], f[(i + 1) % f.size()])];
    std::vector<char> seam(n, 0);
    for (const auto& [k, c] : edgeCount) {
        if (c != 1) continue;
        for (int v : {(int)(k >> 32), (int)(k & 0xffffffff)})
            if (std::abs(out.vertices[v][axis]) <= seamTolerance) seam[v] = 1;
    }

    std::vector<int> mirror(n);
    for (size_t v = 0; v < n; ++v) {
        if (seam[v]) {
            out.vertices[v][axis] = 0.0;
            mirror[v] = (int)v;
        } else {
            auto p = half.vertices[v];
            p[axis] = -p[axis];
            mirror[v] = (int)out.vertices.size();
            out.vertices.push_back(p);
        }
    }
    for (const auto& f : half.faces) {
        std::vector<int> g(f.rbegin(), f.rend());
        for (int& v : g) v = mirror[v];
        out.faces.push_back(std::move(g));
    }
    return out;
}

}  // namespace remersher
