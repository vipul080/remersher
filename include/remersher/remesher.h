#pragma once

#include <functional>
#include <string>

#include "remersher/mesh.h"

namespace remersher {

enum class Resample {
    Auto,    // resample when the input has many sliver triangles
    Always,  // always resample to an isotropic triangulation first
    Never,
};

struct Settings {
    // Desired number of quads in the output.
    int targetQuadCount = 5000;
    // Number of extra solves used to steer the result towards targetQuadCount (0 = single solve).
    int countCalibrationPasses = 2;
    // Stop calibrating once the result is within this fraction of the target.
    double countTolerance = 0.03;
    // Let quads get smaller in high-curvature regions and larger in flat regions.
    bool adaptiveSize = true;
    // Align edge loops to sharp creases of the input.
    bool detectHardEdges = true;
    // Keep open borders of the input as borders of the output.
    bool preserveBoundary = true;
    // Isotropic resampling of the input before solving. The field solver is unstable on inputs
    // with many slivers (e.g. dense UV spheres); resampling fixes that at some extra cost.
    Resample resample = Resample::Auto;
    // Seed for the randomized parts of the solver; results are deterministic per seed.
    int seed = 0;
    // Per-solve time limit. On POSIX each solve runs in a child process, so a solve that hangs or
    // crashes is killed and retried with another seed. 0 disables isolation.
    double solveTimeoutSeconds = 120.0;
    // Retries allowed after failed solves (hung, crashed, or rejected by the quality check), in
    // addition to calibration passes.
    int maxFailedSolves = 4;
};

struct Report {
    int solves = 0;
    int failedSolves = 0;
    int requestedFaces = 0;  // face budget passed to the final (kept) solve
    double seconds = 0.0;
};

using LogFn = std::function<void(const std::string&)>;

// Remeshes a polygon mesh into a quad-dominant mesh. Ngons in the input are fan-triangulated.
// Throws std::runtime_error if the input is invalid or the solver produces no faces.
Mesh remesh(const Mesh& input, const Settings& settings, Report* report = nullptr,
            const LogFn& log = nullptr);

}  // namespace remersher
