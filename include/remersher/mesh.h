#pragma once

#include <array>
#include <string>
#include <vector>

namespace remersher {

// Polygon mesh with 0-based vertex indices. Faces may be triangles, quads or ngons.
struct Mesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::vector<int>> faces;
    // Optional per-vertex density paint in [0, 1] (0.5 = neutral), e.g. from vertex colors.
    // Empty when the input has none.
    std::vector<float> density;

    size_t countFacesWithSides(size_t n) const;
};

// Reads positions and faces from a Wavefront OBJ file. Vertex colors written as `v x y z r g b`
// (a common extension, used by Blender's exporter) become `density` as their luminance.
// Texture/normal indices are ignored, negative (relative) indices are supported. Throws
// std::runtime_error on failure.
Mesh readObj(const std::string& path);

// Writes positions and faces to a Wavefront OBJ file. Throws std::runtime_error on failure.
void writeObj(const std::string& path, const Mesh& mesh);

}  // namespace remersher
