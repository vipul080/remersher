#include "curvature.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace remersher {

using qflow::Hierarchy;
using qflow::Vector3d;

namespace {

struct CurvatureField {
    std::vector<std::vector<int>> nbrs;
    std::vector<Eigen::Matrix3d> tensor;  // smoothed shape operator as a 3D tangent tensor
    std::vector<char> valid;
};

CurvatureField estimateCurvature(const Hierarchy& h) {
    CurvatureField cf;
    const auto& V = h.mV[0];
    const auto& N = h.mN[0];
    const auto& F = h.mF;
    const int n = (int)V.cols();
    auto& nbrs = cf.nbrs;
    auto& tensor = cf.tensor;
    auto& valid = cf.valid;
    nbrs.assign(n, {});
    for (int f = 0; f < F.cols(); ++f)
        for (int i = 0; i < 3; ++i) {
            int a = F(i, f), b = F((i + 1) % 3, f);
            nbrs[a].push_back(b);
            nbrs[b].push_back(a);
        }
    for (auto& nb : nbrs) {
        std::sort(nb.begin(), nb.end());
        nb.erase(std::unique(nb.begin(), nb.end()), nb.end());
    }

    // Per-vertex shape operator, fitted by least squares to the normal curvatures along the edges
    // to its neighbours (kn = 2 n.(vj - vi) / |vj - vi|^2), stored as a 3D tangent tensor so
    // neighbouring estimates can be averaged.
    tensor.assign(n, Eigen::Matrix3d::Zero());
    valid.assign(n, 0);
    for (int v = 0; v < n; ++v) {
        const Vector3d nv = N.col(v);
        if (nbrs[v].size() < 3 || nv.squaredNorm() < 0.5) continue;
        Vector3d u = std::abs(nv.x()) < 0.9 ? Vector3d(1, 0, 0) : Vector3d(0, 1, 0);
        u = (u - nv * nv.dot(u)).normalized();
        const Vector3d w = nv.cross(u);
        Eigen::Matrix3d AtA = Eigen::Matrix3d::Zero();
        Eigen::Vector3d Atb = Eigen::Vector3d::Zero();
        for (int j : nbrs[v]) {
            const Vector3d e = V.col(j) - V.col(v);
            const double len2 = e.squaredNorm();
            Vector3d t = e - nv * nv.dot(e);
            const double tl = t.norm();
            if (len2 <= 0 || tl <= 0) continue;
            t /= tl;
            const double kn = 2.0 * nv.dot(e) / len2;
            const double tx = t.dot(u), ty = t.dot(w);
            const Eigen::Vector3d row(tx * tx, 2 * tx * ty, ty * ty);
            AtA += row * row.transpose();
            Atb += row * kn;
        }
        AtA += Eigen::Matrix3d::Identity() * 1e-9 * AtA.trace();
        const Eigen::Vector3d k = AtA.ldlt().solve(Atb);
        if (!k.allFinite()) continue;
        tensor[v] = k[0] * u * u.transpose() + k[1] * (u * w.transpose() + w * u.transpose()) +
                    k[2] * w * w.transpose();
        valid[v] = 1;
    }

    // Neighbourhood averaging over roughly one target quad, to suppress tessellation noise.
    for (int round = 0; round < 4; ++round) {
        std::vector<Eigen::Matrix3d> next = tensor;
        for (int v = 0; v < n; ++v) {
            if (!valid[v]) continue;
            int count = 1;
            for (int j : nbrs[v])
                if (valid[j]) next[v] += tensor[j], ++count;
            next[v] /= count;
        }
        tensor.swap(next);
    }

    return cf;
}

}  // namespace

int addCurvatureConstraints(Hierarchy& h, double strength) {
    if (strength <= 0 || h.mV.empty()) return 0;
    const auto& N = h.mN[0];
    const int n = (int)h.mV[0].cols();
    CurvatureField cf = estimateCurvature(h);
    const auto& nbrs = cf.nbrs;
    const auto& tensor = cf.tensor;
    const auto& valid = cf.valid;

    // Principal direction and weight per vertex.
    std::vector<Vector3d> dir(n, Vector3d::Zero());
    std::vector<double> weight(n, 0.0);
    for (int v = 0; v < n; ++v) {
        if (!valid[v]) continue;
        const Vector3d nv = N.col(v);
        Vector3d u = std::abs(nv.x()) < 0.9 ? Vector3d(1, 0, 0) : Vector3d(0, 1, 0);
        u = (u - nv * nv.dot(u)).normalized();
        const Vector3d w = nv.cross(u);
        const double a = u.dot(tensor[v] * u), b = u.dot(tensor[v] * w), c = w.dot(tensor[v] * w);
        const double mean = 0.5 * (a + c), diff = std::sqrt(0.25 * (a - c) * (a - c) + b * b);
        const double k1 = mean + diff, k2 = mean - diff;
        const double anisotropy = (k1 - k2) / (std::abs(k1) + std::abs(k2) + 1e-12);
        // Curvature difference over one target quad: how much the direction choice matters here.
        const double significance = (k1 - k2) * h.mScale;
        weight[v] = std::clamp((anisotropy - 0.4) / 0.4, 0.0, 1.0) *
                    std::clamp(significance / 0.15, 0.0, 1.0);
        const double theta = 0.5 * std::atan2(2 * b, a - c);
        dir[v] = std::cos(theta) * u + std::sin(theta) * w;
    }

    // Only constrain where the direction is coherent with the neighbours (as a 4-RoSy: angles
    // compared modulo 90 degrees). Where it swirls, following it would just add singularities.
    int constrained = 0;
    std::vector<double> final(n, 0.0);
    for (int v = 0; v < n; ++v) {
        if (weight[v] <= 0 || h.mCQw[0][v] != 0) continue;  // keep border constraints
        const Vector3d nv = N.col(v);
        const Vector3d perp = nv.cross(dir[v]);
        double coherence = 0;
        int count = 0;
        for (int j : nbrs[v]) {
            if (weight[j] <= 0) continue;
            const double x = dir[j].dot(dir[v]), y = dir[j].dot(perp);
            coherence += std::cos(4 * std::atan2(y, x));
            ++count;
        }
        if (count < 3) continue;
        coherence /= count;
        final[v] = weight[v] * std::clamp((coherence - 0.7) / 0.25, 0.0, 1.0);
    }

    // Erode once (each weight becomes the minimum over its 1-ring), so only regions coherent
    // over a few quads keep a constraint: tubes and tori do, while small bumps whose directions
    // change within a quad or two do not (aligning to those only adds singularities).
    for (int round = 0; round < 1; ++round) {
        std::vector<double> next = final;
        for (int v = 0; v < n; ++v)
            for (int j : nbrs[v]) next[v] = std::min(next[v], final[j]);
        final.swap(next);
    }
    for (int v = 0; v < n; ++v) {
        const double w = strength * final[v];
        if (w < 0.01 || h.mCQw[0][v] != 0) continue;
        h.mCQ[0].col(v) = dir[v];
        h.mCQw[0][v] = w;
        ++constrained;
    }
    return constrained;
}

void applyCurvatureSizing(Hierarchy& h, double adaptivity) {
    if (adaptivity <= 0 || h.mV.empty() || h.mS.empty()) return;
    const auto& N = h.mN[0];
    const int n = (int)h.mV[0].cols();
    CurvatureField cf = estimateCurvature(h);

    // Curvature magnitude sqrt(k1^2 + k2^2) per vertex.
    std::vector<double> kappa(n, 0.0);
    double meanKappa = 0;
    int counted = 0;
    for (int v = 0; v < n; ++v) {
        if (!cf.valid[v]) continue;
        const Vector3d nv = N.col(v);
        Vector3d u = std::abs(nv.x()) < 0.9 ? Vector3d(1, 0, 0) : Vector3d(0, 1, 0);
        u = (u - nv * nv.dot(u)).normalized();
        const Vector3d w = nv.cross(u);
        const double a = u.dot(cf.tensor[v] * u), b = u.dot(cf.tensor[v] * w), c = w.dot(cf.tensor[v] * w);
        kappa[v] = std::sqrt(a * a + 2 * b * b + c * c);
        meanKappa += kappa[v];
        ++counted;
    }
    if (!counted) return;
    meanKappa /= counted;
    // Flat (or nearly flat) input: nothing to adapt to, and the ratios below would be 0/0.
    if (!(meanKappa * h.mScale > 1e-6)) return;
    // Relative quad size (kappa_mean / kappa)^gamma: at slider 100 the edge length is inversely
    // proportional to curvature, at 75 to its square root. The quarter-mean offset keeps flat
    // regions finite; sizes are clamped to 1/4..4 of the base size.
    // Up to 50 (the default) the solver's own mild adaptivity is used as is; 50..100 blends in
    // curvature-driven sizing. Applying it at the default distorted quads across the benchmark.
    const double gamma = std::clamp((adaptivity - 50.0) / 50.0, 0.0, 1.0);
    if (gamma <= 0) return;
    std::vector<double> logSize(n, 0.0);
    for (int v = 0; v < n; ++v) {
        if (!cf.valid[v]) continue;
        const double ratio = (meanKappa * 1.25) / (kappa[v] + 0.25 * meanKappa);
        logSize[v] = std::clamp(gamma * std::log(ratio), std::log(0.25), std::log(4.0));
    }
    // Smooth in log space so the size changes gradually (quad meshes need room to change density).
    for (int round = 0; round < 8; ++round) {
        std::vector<double> next = logSize;
        for (int v = 0; v < n; ++v) {
            double sum = logSize[v];
            for (int j : cf.nbrs[v]) sum += logSize[j];
            next[v] = sum / (1 + cf.nbrs[v].size());
        }
        logSize.swap(next);
    }
    std::vector<double> size(n);
    for (int v = 0; v < n; ++v) size[v] = std::exp(logSize[v]);
    // Keep the mean at 1 so the face budget still means the same thing.
    double mean = 0;
    for (double x : size) mean += x;
    mean /= std::max(1, n);
    auto& S = h.mS[0];
    for (int v = 0; v < n; ++v) {
        S(0, v) *= size[v] / mean;
        S(1, v) *= size[v] / mean;
    }
    for (size_t l = 0; l + 1 < h.mS.size(); ++l) {
        const auto& toUpper = h.mToUpper[l];
        for (int i = 0; i < toUpper.cols(); ++i) {
            const int u0 = toUpper(0, i), u1 = toUpper(1, i);
            if (u1 != -1) h.mS[l + 1].col(i) = 0.5 * (h.mS[l].col(u0) + h.mS[l].col(u1));
            else h.mS[l + 1].col(i) = h.mS[l].col(u0);
        }
    }
}

}  // namespace remersher
