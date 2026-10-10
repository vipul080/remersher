#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <string>

#include "remersher/mesh.h"
#include "remersher/remesher.h"

namespace {

const char* kUsage =
    "usage: remersher -i <input.obj> -o <output.obj> [options]\n"
    "\n"
    "options:\n"
    "  -t, --target <n>         target quad count (default 5000)\n"
    "      --passes <n>         extra solves to hit the target count (default 6)\n"
    "  -j, --jobs <n>           calibration solves run in parallel (default 3)\n"
    "      --tolerance <f>      acceptable count error as a fraction (default 0.03)\n"
    "      --no-adaptive        uniform quad size instead of curvature-adaptive\n"
    "      --no-hard-edges      do not align edge loops to sharp creases\n"
    "      --curvature <f>      edge-loop alignment to curvature, 0..1 (default 0.5)\n"
    "      --no-boundary        do not constrain open borders\n"
    "      --resample <mode>    auto | always | never: isotropic resampling of the input\n"
    "                           before solving (default auto: only for sliver-heavy input)\n"
    "      --seed <n>           solver seed (default 0)\n"
    "      --timeout <s>        per-solve time limit in seconds, 0 = none (default 120)\n"
    "  -q, --quiet              only print errors\n"
    "  -h, --help               show this help\n";

[[noreturn]] void usageError(const std::string& msg) {
    std::fprintf(stderr, "remersher: %s\n\n%s", msg.c_str(), kUsage);
    std::exit(2);
}

}  // namespace

int main(int argc, char** argv) {
    std::string input, output;
    remersher::Settings settings;
    bool quiet = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= argc) usageError("missing value for " + arg);
            return argv[++i];
        };
        if (arg == "-i" || arg == "--input") input = value();
        else if (arg == "-o" || arg == "--output") output = value();
        else if (arg == "-t" || arg == "--target") settings.targetQuadCount = std::atoi(value().c_str());
        else if (arg == "--passes") settings.countCalibrationPasses = std::atoi(value().c_str());
        else if (arg == "-j" || arg == "--jobs") settings.maxParallelSolves = std::atoi(value().c_str());
        else if (arg == "--tolerance") settings.countTolerance = std::atof(value().c_str());
        else if (arg == "--no-adaptive") settings.adaptiveSize = false;
        else if (arg == "--no-hard-edges") settings.detectHardEdges = false;
        else if (arg == "--curvature") settings.curvatureAlignment = std::atof(value().c_str());
        else if (arg == "--no-boundary") settings.preserveBoundary = false;
        else if (arg == "--resample") {
            std::string mode = value();
            if (mode == "auto") settings.resample = remersher::Resample::Auto;
            else if (mode == "always") settings.resample = remersher::Resample::Always;
            else if (mode == "never") settings.resample = remersher::Resample::Never;
            else usageError("--resample must be auto, always or never");
        }
        else if (arg == "--seed") settings.seed = std::atoi(value().c_str());
        else if (arg == "--timeout") settings.solveTimeoutSeconds = std::atof(value().c_str());
        else if (arg == "-q" || arg == "--quiet") quiet = true;
        else if (arg == "-h" || arg == "--help") { std::fputs(kUsage, stdout); return 0; }
        else usageError("unknown option " + arg);
    }
    if (input.empty() || output.empty()) usageError("both -i and -o are required");
    if (settings.targetQuadCount <= 0) usageError("--target must be positive");

    try {
        remersher::Mesh mesh = remersher::readObj(input);
        remersher::Report report;
        remersher::LogFn log;
        const auto start = std::chrono::steady_clock::now();
        if (!quiet)
            log = [start](const std::string& line) {
                double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                std::fprintf(stderr, "[%6.2fs] %s\n", t, line.c_str());
            };
        remersher::Mesh result = remersher::remesh(mesh, settings, &report, log);
        remersher::writeObj(output, result);
        if (!quiet)
            std::fprintf(stderr, "wrote %s: %zu faces (%zu quads, %zu tris) in %.2fs, %d solve(s), %d failed\n",
                         output.c_str(), result.faces.size(), result.countFacesWithSides(4),
                         result.countFacesWithSides(3), report.seconds, report.solves, report.failedSolves);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "remersher: %s\n", e.what());
        return 1;
    }
    return 0;
}
