#include "remersher/remesher.h"

#include "curvature.h"
#include "crease.h"
#include "density.h"
#include "feature_lines.h"
#include "features.h"
#include "mesh_cleanup.h"
#include "quality.h"
#include "resample.h"
#include "symmetry.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <cstring>
#include <memory>
#include <stdexcept>

#include "config.hpp"
#include "field-math.hpp"
#include "optimizer.hpp"
#include "parametrizer.hpp"

#if defined(__unix__) || defined(__APPLE__)
#define REMERSHER_ISOLATE_SOLVES 1
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <poll.h>
#include <sys/wait.h>
#include <unistd.h>
#else
#define REMERSHER_ISOLATE_SOLVES 0
#endif

namespace remersher {
namespace {

using Clock = std::chrono::steady_clock;

// Triangle mesh in QuadriFlow's layout (3 x N columns).
struct TriangleInput {
    qflow::MatrixXd V;
    qflow::MatrixXi F;
    const DensitySource* density = nullptr;  // optional density paint
    const FeatureLines* features = nullptr;  // optional forced feature lines (input coordinates)
};

TriangleInput toSolverInput(const TriangleMesh& mesh) {
    TriangleInput out;
    out.V.resize(3, mesh.vertices.size());
    for (size_t i = 0; i < mesh.vertices.size(); ++i)
        for (int j = 0; j < 3; ++j) out.V(j, i) = mesh.vertices[i][j];
    out.F.resize(3, mesh.triangles.size());
    for (size_t i = 0; i < mesh.triangles.size(); ++i)
        for (int j = 0; j < 3; ++j) out.F(j, i) = mesh.triangles[i][j];
    return out;
}

// Converts QuadriFlow's compact quad output to a Mesh: collapses repeated corners (quads that
// degenerated to triangles), drops faces with fewer than 3 distinct corners and unused vertices.
Mesh fromParametrizer(const qflow::Parametrizer& field) {
    Mesh out;
    std::vector<int> remap(field.O_compact.size(), -1);
    for (const auto& q : field.F_compact) {
        std::vector<int> face;
        for (int j = 0; j < 4; ++j) {
            int v = q[j];
            if (v < 0 || v >= (int)remap.size()) { face.clear(); break; }
            if (std::find(face.begin(), face.end(), v) == face.end()) face.push_back(v);
        }
        if (face.size() < 3) continue;
        for (int& v : face) {
            if (remap[v] < 0) {
                remap[v] = (int)out.vertices.size();
                auto p = field.O_compact[v] * field.normalize_scale + field.normalize_offset;
                out.vertices.push_back({p[0], p[1], p[2]});
            }
            v = remap[v];
        }
        out.faces.push_back(std::move(face));
    }
    return out;
}

Mesh solveOnce(const TriangleInput& input, const Settings& s, int faceBudget) {
    // Parametrizer is large and keeps per-solve state, so every solve gets a fresh one.
    auto field = std::make_unique<qflow::Parametrizer>();
    field->V = input.V;
    field->F = input.F;
    field->NormalizeMesh();
    field->flag_preserve_sharp = s.detectHardEdges ? 1 : 0;
    field->sharp_angle_degrees = s.hardEdgeAngle;
    field->flag_preserve_boundary = s.preserveBoundary ? 1 : 0;
    field->flag_adaptive_scale = s.adaptivity > 0 ? 1 : 0;
    field->hierarchy.rng_seed = s.seed;

    field->Initialize(faceBudget);

    qflow::Hierarchy& h = field->hierarchy;
    // Forced feature lines: half-edges of the solver mesh lying on them, in input coordinates.
    std::vector<char> forcedEdge;
    if (input.features) {
        forcedEdge.assign(3 * h.mF.cols(), 0);
        auto toInput = [&](int v) {
            qflow::Vector3d p = h.mV[0].col(v) * field->normalize_scale + field->normalize_offset;
            return std::array<double, 3>{p[0], p[1], p[2]};
        };
        for (uint32_t i = 0; i < 3 * h.mF.cols(); ++i) {
            if (!input.features->covers(toInput(h.mF(i % 3, i / 3)), toInput(h.mF((i + 1) % 3, i / 3)))) continue;
            forcedEdge[i] = 1;
            field->sharp_edges[i] = 1;  // the integer stage keeps these as edges
            if (h.mE2E[i] >= 0) forcedEdge[h.mE2E[i]] = field->sharp_edges[h.mE2E[i]] = 1;
        }
    }
    h.clearConstraints();
    if (field->flag_preserve_boundary) {
        // Constrain the orientation and position fields to follow open borders.
        for (uint32_t i = 0; i < 3 * h.mF.cols(); ++i) {
            if (h.mE2E[i] != -1) continue;
            uint32_t i0 = h.mF(i % 3, i / 3);
            uint32_t i1 = h.mF((i + 1) % 3, i / 3);
            qflow::Vector3d p0 = h.mV[0].col(i0), p1 = h.mV[0].col(i1);
            qflow::Vector3d edge = p1 - p0;
            if (edge.squaredNorm() <= 0) continue;
            edge.normalize();
            h.mCO[0].col(i0) = p0;
            h.mCO[0].col(i1) = p1;
            h.mCQ[0].col(i0) = h.mCQ[0].col(i1) = edge;
            h.mCQw[0][i0] = h.mCQw[0][i1] = h.mCOw[0][i0] = h.mCOw[0][i1] = 1.0;
        }
    }
    if (field->flag_preserve_sharp || !forcedEdge.empty()) {
        // Constrain the fields to run along hard edges, as QuadriFlow does for borders. Without
        // this the orientation field ignores creases and only the later integer stage tries to
        // pull edge loops onto them, which leaves misplaced corner singularities.
        const double cosHard = field->flag_preserve_sharp ? std::cos(s.hardEdgeAngle * M_PI / 180.0) : 2.0;
        std::vector<std::vector<qflow::Vector3d>> dirs(h.mV[0].cols());
        auto faceNormal = [&](uint32_t e) {
            const uint32_t f = e / 3;
            qflow::Vector3d a = h.mV[0].col(h.mF(0, f)), b = h.mV[0].col(h.mF(1, f)), c = h.mV[0].col(h.mF(2, f));
            return (b - a).cross(c - a).normalized();
        };
        for (uint32_t i = 0; i < 3 * h.mF.cols(); ++i) {
            const int j = h.mE2E[i];
            if (j < 0 || (uint32_t)j < i) continue;  // borders handled above; each pair once
            const bool forced = !forcedEdge.empty() && forcedEdge[i];
            if (!forced && faceNormal(i).dot(faceNormal(j)) >= cosHard) continue;
            const uint32_t i0 = h.mF(i % 3, i / 3), i1 = h.mF((i + 1) % 3, i / 3);
            const qflow::Vector3d edge = h.mV[0].col(i1) - h.mV[0].col(i0);
            if (edge.squaredNorm() <= 0) continue;
            dirs[i0].push_back(edge.normalized());
            dirs[i1].push_back(edge.normalized());
        }
        for (int v = 0; v < (int)dirs.size(); ++v) {
            if (dirs[v].empty() || h.mCQw[0][v] != 0) continue;
            const qflow::Vector3d nv = h.mN[0].col(v);
            qflow::Vector3d d0 = dirs[v][0] - nv * nv.dot(dirs[v][0]);
            if (d0.squaredNorm() < 1e-12) continue;
            d0.normalize();
            // Only where the creases through this vertex agree as a 4-RoSy (meet at multiples of
            // 90 degrees). At other corners (gear teeth, cube corners) the field needs a
            // singularity, so it is left free to place one.
            const qflow::Vector3d perp = nv.cross(d0);
            bool consistent = true;
            for (size_t k = 1; k < dirs[v].size() && consistent; ++k) {
                const double x = dirs[v][k].dot(d0), y = dirs[v][k].dot(perp);
                consistent = std::cos(4 * std::atan2(y, x)) > 0.5;  // within ~15 degrees
            }
            if (!consistent) continue;
            h.mCQ[0].col(v) = d0;
            h.mCO[0].col(v) = h.mV[0].col(v);
            h.mCQw[0][v] = h.mCOw[0][v] = 1.0;
        }
    }
    addCurvatureConstraints(h, s.curvatureAlignment);
    h.propagateConstraints();

    qflow::Optimizer::optimize_orientations(field->hierarchy);
    field->ComputeOrientationSingularities();
    if (field->flag_adaptive_scale) field->EstimateSlope();
    qflow::Optimizer::optimize_scale(field->hierarchy, field->rho, field->flag_adaptive_scale);
    // QuadriFlow's own "adaptive" field only varies by about +-25% (it mostly follows the field's
    // slope); add a curvature-driven size on top so the slider behaves like users expect.
    applyCurvatureSizing(field->hierarchy, s.adaptivity);
    if (input.density)
        applyDensityPaint(field->hierarchy, *input.density, field->normalize_scale, field->normalize_offset);
    field->flag_adaptive_scale = 1;  // the position solve always uses the scale field
    qflow::Optimizer::optimize_positions(field->hierarchy, field->flag_adaptive_scale);
    field->ComputePositionSingularities();
    // Upstream calls this without the sizing field, which made "adaptive" output uniform.
    field->ComputeIndexMap(s.adaptivity > 50 || input.density ? 1 : 0);

    return fromParametrizer(*field);
}

#if REMERSHER_ISOLATE_SOLVES
// Mesh wire format between the solve child and the parent:
// [u64 nv][u64 nf][nv * 3 doubles][nf * (u32 size, size * i32 indices)]
void appendBytes(std::string& buf, const void* p, size_t n) { buf.append((const char*)p, n); }

std::string serialize(const Mesh& m) {
    std::string buf;
    uint64_t nv = m.vertices.size(), nf = m.faces.size();
    appendBytes(buf, &nv, 8);
    appendBytes(buf, &nf, 8);
    for (const auto& v : m.vertices) appendBytes(buf, v.data(), 24);
    for (const auto& f : m.faces) {
        uint32_t n = (uint32_t)f.size();
        appendBytes(buf, &n, 4);
        appendBytes(buf, f.data(), 4 * n);
    }
    return buf;
}

bool deserialize(const std::string& buf, Mesh& m) {
    size_t pos = 0;
    auto take = [&](void* dst, size_t n) {
        if (pos + n > buf.size()) return false;
        std::memcpy(dst, buf.data() + pos, n);
        pos += n;
        return true;
    };
    uint64_t nv = 0, nf = 0;
    if (!take(&nv, 8) || !take(&nf, 8) || nv > buf.size() || nf > buf.size()) return false;
    m.vertices.resize(nv);
    for (auto& v : m.vertices)
        if (!take(v.data(), 24)) return false;
    m.faces.resize(nf);
    for (auto& f : m.faces) {
        uint32_t n = 0;
        if (!take(&n, 4) || n > nv) return false;
        f.resize(n);
        if (!take(f.data(), 4 * n)) return false;
        for (int idx : f)
            if (idx < 0 || (uint64_t)idx >= nv) return false;
    }
    return pos == buf.size();
}

#endif

struct SolveJob {
    Settings settings;
    int budget = 0;
};

struct SolveOutcome {
    bool ok = false;
    Mesh mesh;
    std::string why;
};

bool isUsable(const Mesh& m, std::string& why) {
    for (const auto& v : m.vertices) {
        if (!std::isfinite(v[0]) || !std::isfinite(v[1]) || !std::isfinite(v[2])) {
            why = "non-finite vertex positions";
            return false;
        }
    }
    return true;
}

#if REMERSHER_ISOLATE_SOLVES
// Runs each job's solveOnce in its own forked child, all concurrently, so hangs and crashes inside
// the solver cannot take down the caller and calibration candidates use several cores.
std::vector<SolveOutcome> solveIsolated(const TriangleInput& input, const std::vector<SolveJob>& jobs,
                                        double timeout) {
    struct Child {
        pid_t pid = -1;
        int fd = -1;
        std::string buf;
        bool done = false;
    };
    std::vector<SolveOutcome> out(jobs.size());
    std::vector<Child> kids(jobs.size());
    std::fflush(nullptr);
    for (size_t i = 0; i < jobs.size(); ++i) {
        int fds[2];
        if (pipe(fds) != 0) { out[i].why = "pipe failed"; kids[i].done = true; continue; }
        pid_t pid = fork();
        if (pid < 0) {
            close(fds[0]);
            close(fds[1]);
            out[i].why = "fork failed";
            kids[i].done = true;
            continue;
        }
        if (pid == 0) {
            close(fds[0]);
            int rc = 0;
            try {
                std::string buf = serialize(solveOnce(input, jobs[i].settings, jobs[i].budget));
                for (size_t off = 0; off < buf.size();) {
                    ssize_t n = write(fds[1], buf.data() + off, buf.size() - off);
                    if (n <= 0) { rc = 3; break; }
                    off += (size_t)n;
                }
            } catch (...) {
                rc = 2;
            }
            close(fds[1]);
            _exit(rc);
        }
        close(fds[1]);
        kids[i].pid = pid;
        kids[i].fd = fds[0];
    }

    char chunk[1 << 16];
    const auto deadline = Clock::now() + std::chrono::duration<double>(timeout);
    for (;;) {
        std::vector<pollfd> pfds;
        std::vector<size_t> owner;
        for (size_t i = 0; i < kids.size(); ++i)
            if (!kids[i].done) pfds.push_back({kids[i].fd, POLLIN, 0}), owner.push_back(i);
        if (pfds.empty()) break;
        double left = std::chrono::duration<double>(deadline - Clock::now()).count();
        if (left <= 0) break;
        int ready = poll(pfds.data(), pfds.size(), (int)std::min(left * 1000.0 + 1, 1000.0));
        if (ready < 0 && errno != EINTR) break;
        if (ready <= 0) continue;
        for (size_t k = 0; k < pfds.size(); ++k) {
            if (!(pfds[k].revents & (POLLIN | POLLHUP | POLLERR))) continue;
            Child& c = kids[owner[k]];
            ssize_t n = read(c.fd, chunk, sizeof chunk);
            if (n > 0) c.buf.append(chunk, (size_t)n);
            else if (n == 0 || errno != EINTR) c.done = true;  // EOF or error
        }
    }

    for (size_t i = 0; i < kids.size(); ++i) {
        Child& c = kids[i];
        if (c.pid < 0) continue;
        const bool timedOut = !c.done;
        close(c.fd);
        if (timedOut) kill(c.pid, SIGKILL);
        int status = 0;
        while (waitpid(c.pid, &status, 0) < 0 && errno == EINTR) {}
        SolveOutcome& o = out[i];
        if (timedOut) o.why = "timed out";
        else if (WIFSIGNALED(status)) o.why = std::string("crashed (") + strsignal(WTERMSIG(status)) + ")";
        else if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) o.why = "solver error";
        else if (!deserialize(c.buf, o.mesh)) o.why = "bad result";
        else o.ok = isUsable(o.mesh, o.why);
    }
    return out;
}
#endif

std::vector<SolveOutcome> solveMany(const TriangleInput& input, const std::vector<SolveJob>& jobs) {
#if REMERSHER_ISOLATE_SOLVES
    if (!jobs.empty() && jobs[0].settings.solveTimeoutSeconds > 0)
        return solveIsolated(input, jobs, jobs[0].settings.solveTimeoutSeconds);
#endif
    std::vector<SolveOutcome> out(jobs.size());
    for (size_t i = 0; i < jobs.size(); ++i) {
        try {
            out[i].mesh = solveOnce(input, jobs[i].settings, jobs[i].budget);
            out[i].ok = isUsable(out[i].mesh, out[i].why);
        } catch (const std::exception& e) {
            out[i].why = e.what();
        }
    }
    return out;
}

}  // namespace

Mesh remesh(const Mesh& input, const Settings& settings, Report* report, const LogFn& log) {
    if (settings.targetQuadCount <= 0) throw std::runtime_error("targetQuadCount must be positive");
    const auto start = Clock::now();

    // Forced feature lines: explicit `l` edges plus, if requested, borders between materials.
    std::vector<std::array<FeatureLines::Vec3, 2>> featureSegments;
    for (const auto& e : input.featureEdges)
        if (e[0] != e[1]) featureSegments.push_back({input.vertices[e[0]], input.vertices[e[1]]});
    if (settings.useMaterials && input.faceMaterial.size() == input.faces.size()) {
        std::map<std::pair<int, int>, int> edgeMaterial;
        for (size_t fi = 0; fi < input.faces.size(); ++fi) {
            const auto& f = input.faces[fi];
            for (size_t i = 0; i < f.size(); ++i) {
                const std::pair<int, int> k = std::minmax(f[i], f[(i + 1) % f.size()]);
                auto [it, inserted] = edgeMaterial.emplace(k, input.faceMaterial[fi]);
                if (!inserted && it->second != input.faceMaterial[fi] && it->second != -2) {
                    featureSegments.push_back({input.vertices[k.first], input.vertices[k.second]});
                    it->second = -2;  // record each border edge once
                }
            }
        }
    }
    double diag = 0;
    if (!input.vertices.empty()) {
        std::array<double, 3> lo = input.vertices[0], hi = lo;
        for (const auto& v : input.vertices)
            for (int j = 0; j < 3; ++j) lo[j] = std::min(lo[j], v[j]), hi[j] = std::max(hi[j], v[j]);
        diag = std::sqrt((hi[0] - lo[0]) * (hi[0] - lo[0]) + (hi[1] - lo[1]) * (hi[1] - lo[1]) + (hi[2] - lo[2]) * (hi[2] - lo[2]));
    }
    const FeatureLines forced(featureSegments, 1e-5 * diag);
    const FeatureLines* forcedPtr = forced.empty() ? nullptr : &forced;
    if (log && forcedPtr) log("feature lines: " + std::to_string(forced.segments().size()) + " edges");

    if (settings.symmetryAxes) {
        // A mesh that does not cross a symmetry plane through the origin (e.g. modelled away from
        // its origin) is mirrored about its bounding-box centre instead; mirroring about the
        // origin would add a detached copy.
        std::array<double, 3> shift{0, 0, 0};
        bool shifted = false;
        for (int a = 0; a < 3; ++a) {
            if (!(settings.symmetryAxes & (1 << a)) || input.vertices.empty()) continue;
            double lo = INFINITY, hi = -INFINITY;
            for (const auto& v : input.vertices) lo = std::min(lo, v[a]), hi = std::max(hi, v[a]);
            const double extent = hi - lo;
            if (lo < -0.05 * extent && hi > 0.05 * extent) continue;
            shift[a] = 0.5 * (lo + hi);
            shifted = true;
        }
        if (shifted) {
            if (log) log("symmetry: mesh does not cross the origin plane; mirroring about its centre");
            Mesh moved = input;
            for (auto& v : moved.vertices)
                for (int a = 0; a < 3; ++a) v[a] -= shift[a];
            Mesh result = remesh(moved, settings, report, log);
            for (auto& v : result.vertices)
                for (int a = 0; a < 3; ++a) v[a] += shift[a];
            return result;
        }
        TriangleMesh half = cleanupForRemeshing(input, 1e-6);
        int axes = 0;
        for (int a = 0; a < 3; ++a)
            if (settings.symmetryAxes & (1 << a)) half = clipToHalfSpace(half, a), ++axes;
        if (half.triangles.empty()) throw std::runtime_error("no geometry on the positive side of the symmetry plane");
        Mesh halfMesh;
        halfMesh.vertices = half.vertices;
        if (settings.useVertexColor && input.density.size() == input.vertices.size()) {
            const DensitySource paint(input.vertices, input.density);
            for (const auto& p : halfMesh.vertices) halfMesh.density.push_back(paint.sample(p));
        }
        // Feature lines on the kept side (cut at the planes) travel as extra `l` edges.
        for (const auto& seg : featureSegments) {
            std::array<double, 3> a = seg[0], b = seg[1];
            bool keep = true;
            for (int ax = 0; ax < 3 && keep; ++ax) {
                if (!(settings.symmetryAxes & (1 << ax))) continue;
                if (a[ax] < 0 && b[ax] < 0) keep = false;
                else if (a[ax] < 0 || b[ax] < 0) {
                    const double t = a[ax] / (a[ax] - b[ax]);
                    std::array<double, 3> c;
                    for (int j = 0; j < 3; ++j) c[j] = a[j] + t * (b[j] - a[j]);
                    c[ax] = 0;
                    (a[ax] < 0 ? a : b) = c;
                }
            }
            if (!keep) continue;
            const int ia = (int)halfMesh.vertices.size();
            halfMesh.vertices.push_back(a);
            halfMesh.vertices.push_back(b);
            if (!halfMesh.density.empty()) halfMesh.density.insert(halfMesh.density.end(), 2, 0.5f);
            halfMesh.featureEdges.push_back({ia, ia + 1});
        }
        for (const auto& t : half.triangles) halfMesh.faces.push_back({t[0], t[1], t[2]});
        Settings sub = settings;
        sub.symmetryAxes = 0;
        sub.useMaterials = false;  // already turned into feature edges above
        sub.preserveBoundary = true;  // the cut must stay a clean border to weld the halves
        sub.targetQuadCount = std::max(1, settings.targetQuadCount >> axes);
        Mesh result = remesh(halfMesh, sub, report, log);
        double edge = 0;
        size_t edges = 0;
        for (const auto& f : result.faces)
            for (size_t i = 0; i < f.size(); ++i, ++edges) {
                const auto &p = result.vertices[f[i]], &q = result.vertices[f[(i + 1) % f.size()]];
                edge += std::sqrt((p[0] - q[0]) * (p[0] - q[0]) + (p[1] - q[1]) * (p[1] - q[1]) + (p[2] - q[2]) * (p[2] - q[2]));
            }
        const double seamTolerance = edges ? 0.3 * edge / edges : 0.0;
        for (int a = 2; a >= 0; --a)
            if (settings.symmetryAxes & (1 << a)) result = mirrorAcross(result, a, seamTolerance);
        repairPolygonMesh(result);
        if (report) report->seconds = std::chrono::duration<double>(Clock::now() - start).count();
        return result;
    }
    CleanupStats cleanup;
    const TriangleMesh clean = cleanupForRemeshing(input, 1e-6, &cleanup);
    // Hard edges are found once, on the cleaned input, and from here on travel as feature lines
    // together with the user's (materials, sharp edges): resampling and the solver's subdivision
    // add flat edges that would make the relative crease test meaningless there.
    std::vector<std::array<FeatureLines::Vec3, 2>> allSegments = featureSegments;
    if (settings.detectHardEdges) {
        // Bevels narrower than about half a quad count as one hard edge at the output's scale.
        double area = 0;
        for (const auto& t : clean.triangles) {
            const auto &a = clean.vertices[t[0]], &b = clean.vertices[t[1]], &c = clean.vertices[t[2]];
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]}, w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            area += 0.5 * std::sqrt(std::pow(u[1] * w[2] - u[2] * w[1], 2) + std::pow(u[2] * w[0] - u[0] * w[2], 2) +
                                    std::pow(u[0] * w[1] - u[1] * w[0], 2));
        }
        const double quadEdge = std::sqrt(area / settings.targetQuadCount);
        auto creases = detectCreases(clean, settings.hardEdgeAngle, 0.5 * quadEdge);
        if (log && !creases.empty()) log("hard edges: " + std::to_string(creases.size()));
        allSegments.insert(allSegments.end(), creases.begin(), creases.end());
    }
    const FeatureLines features(allSegments, 1e-5 * diag);
    const FeatureLines* featuresPtr = features.empty() ? nullptr : &features;
    TriangleMesh solverMesh = clean;
    const double slivers = sliverFraction(clean);
    if (settings.resample == Resample::Always ||
        (settings.resample == Resample::Auto && slivers > 0.05)) {
        double area = 0;
        for (const auto& t : clean.triangles) {
            const auto &a = clean.vertices[t[0]], &b = clean.vertices[t[1]], &c = clean.vertices[t[2]];
            const double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
            const double w[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
            area += 0.5 * std::sqrt(std::pow(u[1] * w[2] - u[2] * w[1], 2) + std::pow(u[2] * w[0] - u[0] * w[2], 2) +
                                    std::pow(u[0] * w[1] - u[1] * w[0], 2));
        }
        // A bit finer than the solver's own working resolution (half the quad size), so
        // adaptive sizing still has room to shrink quads in curved regions.
        const double edge = 0.4 * std::sqrt(area / settings.targetQuadCount);
        ResampleStats rs;
        solverMesh = resampleIsotropic(clean, edge, 0.0, 5, &rs, featuresPtr);
        if (log) {
            char line[200];
            std::snprintf(line, sizeof line,
                          "resample: %.0f%% slivers -> %zu triangles at edge %.4g (%d splits, %d collapses, %d flips)",
                          100 * slivers, solverMesh.triangles.size(), edge, rs.splits, rs.collapses, rs.flips);
            log(line);
        }
    }
    TriangleInput tris = toSolverInput(solverMesh);
    tris.features = featuresPtr;
    std::unique_ptr<DensitySource> paint;
    if (settings.useVertexColor && input.density.size() == input.vertices.size() && !input.vertices.empty()) {
        paint = std::make_unique<DensitySource>(input.vertices, input.density);
        tris.density = paint.get();
        if (log) log("using vertex-color density");
    }
    bool usingResampled = solverMesh.triangles != clean.triangles;
    int consecutiveFailures = 0;
    if (log && (cleanup.weldedVertices || cleanup.droppedTriangles || cleanup.flippedTriangles))
        log("cleanup: welded " + std::to_string(cleanup.weldedVertices) + " vertices, dropped " +
            std::to_string(cleanup.droppedTriangles) + " triangles, re-wound " +
            std::to_string(cleanup.flippedTriangles) + " triangles");
    const double target = settings.targetQuadCount;

    Mesh best;
    double bestError = INFINITY;
    int bestBudget = 0;
    int budget = settings.targetQuadCount;
    int solves = 0, failed = 0;
    Settings attempt = settings;
    attempt.hardEdgeAngle = 180;  // creases come in as feature lines; no angle test in the solver
    // Least-bad result among those rejected by the quality check, returned (with a warning) only
    // if no solve passes: a flawed mesh is more useful to the user than an error.
    Mesh fallback;
    double fallbackBadness = INFINITY;
    int fallbackBudget = 0;

    // Repairs, snaps and quality-checks one solver result. Returns false if it was rejected.
    auto finish = [&](SolveOutcome& o, int budgetUsed) {
        if (!o.ok) return false;
        repairPolygonMesh(o.mesh);
        snapToFeatures(o.mesh, clean, 0.0, settings.relaxIterations, settings.cleanupPoles, featuresPtr);
        repairPolygonMesh(o.mesh);  // snapping can land two vertices on the same point
        // The solver occasionally collapses or folds whole regions without reporting an error;
        // reject those results like a crash.
        const QualityCheck qc = checkQuality(clean, o.mesh);
        if (qc.acceptable(&o.why)) return true;
        const double badness = qc.farOutputFraction + qc.uncoveredFraction + qc.flippedFraction;
        if (!o.mesh.faces.empty() && badness < fallbackBadness) {
            fallbackBadness = badness;
            fallback = o.mesh;
            fallbackBudget = budgetUsed;
        }
        char detail[128];
        std::snprintf(detail, sizeof detail, " (far %.1f%%, holes %.1f%%, folded %.1f%%)",
                      100 * qc.farOutputFraction, 100 * qc.uncoveredFraction, 100 * qc.flippedFraction);
        o.why += detail;
        o.ok = false;
        return false;
    };
    // Every accepted (budget, face count) pair, for bracketing the target during calibration.
    std::vector<std::pair<int, size_t>> tried;
    auto consider = [&](Mesh& m, int budgetUsed) {
        ++solves;
        tried.push_back({budgetUsed, m.faces.size()});
        const double error = std::abs((double)m.faces.size() - target) / target;
        if (log)
            log("solve " + std::to_string(solves) + ": budget " + std::to_string(budgetUsed) + " -> " +
                std::to_string(m.faces.size()) + " faces");
        if (!m.faces.empty() && error < bestError) {
            bestError = error;
            best = std::move(m);
            bestBudget = budgetUsed;
        }
    };

    // 1. A first solve, retried with perturbed seed/budget (and finally without resampling) until
    //    one passes.
    for (;;) {
        std::vector<SolveOutcome> res = solveMany(tris, {{attempt, budget}});
        if (finish(res[0], budget)) {
            consider(res[0].mesh, budget);
            break;
        }
        ++failed;
        ++consecutiveFailures;
        if (log) log("solve rejected at budget " + std::to_string(budget) + ": " + res[0].why);
        if (failed > settings.maxFailedSolves) break;
        if (usingResampled && consecutiveFailures >= 2) {
            // Resampling can occasionally produce input the solver chokes on; fall back to the
            // cleaned original rather than burning the remaining attempts.
            if (log) log("falling back to the original triangulation");
            tris = toSolverInput(clean);
            tris.density = paint.get();
            tris.features = featuresPtr;
            usingResampled = false;
            consecutiveFailures = 0;
            budget = settings.targetQuadCount;
            continue;
        }
        // Failures are specific to one seed/budget combination; perturb both and retry.
        attempt.seed += 1;
        budget = std::max(1, (int)std::lround(budget * 1.03));
    }

    // 2. Count calibration. The face budget only sets the target edge length, so the solver's
    //    count drifts (singularities, adaptivity, hard edges) and is noisy. It scales roughly
    //    linearly with the budget, so solve a few candidates around the corrected budget at once
    //    and keep the closest.
    int passesLeft = std::max(0, settings.countCalibrationPasses);
    while (!best.faces.empty() && bestError > settings.countTolerance && passesLeft > 0) {
        const int n = std::min(passesLeft, std::max(1, settings.maxParallelSolves));
        std::vector<SolveJob> jobs;
        // If earlier solves landed on both sides of the target, search between the closest budgets
        // below and above it: the count can jump with the budget, and a ratio step overshoots.
        int below = -1, above = -1;
        for (const auto& [b, count] : tried) {
            if ((double)count < target && b > below) below = b;
            if ((double)count > target && (above < 0 || b < above)) above = b;
        }
        if (below > 0 && above > below + n) {
            for (int i = 1; i <= n; ++i) jobs.push_back({attempt, below + (above - below) * i / (n + 1)});
        } else {
            const double ratio = std::clamp(target / (double)best.faces.size(), 0.25, 4.0);
            const double estimate = bestBudget * ratio;
            for (int i = 0; i < n; ++i) {
                const double spread = n == 1 ? 1.0 : 1.0 + 0.08 * (i - (n - 1) / 2.0);
                jobs.push_back({attempt, std::max(1, (int)std::lround(estimate * spread))});
            }
        }
        passesLeft -= n;
        std::vector<SolveOutcome> res = solveMany(tris, jobs);
        for (size_t i = 0; i < res.size(); ++i) {
            if (finish(res[i], jobs[i].budget)) {
                consider(res[i].mesh, jobs[i].budget);
            } else {
                ++failed;
                if (log) log("solve rejected at budget " + std::to_string(jobs[i].budget) + ": " + res[i].why);
            }
        }
        attempt.seed += 1;  // a fresh seed if another round follows
    }

    if (best.faces.empty() && !fallback.faces.empty()) {
        if (log) log("warning: no solve passed the quality check; returning the best attempt");
        best = std::move(fallback);
        bestBudget = fallbackBudget;
    }
    if (best.faces.empty())
        throw std::runtime_error(solves ? "solver produced an empty mesh" : "all solves failed");
    if (report) {
        report->solves = solves;
        report->failedSolves = failed;
        report->requestedFaces = bestBudget;
        report->seconds = std::chrono::duration<double>(Clock::now() - start).count();
    }
    return best;
}

}  // namespace remersher
