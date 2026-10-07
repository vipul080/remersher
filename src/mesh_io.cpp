#include "remersher/mesh.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace remersher {

size_t Mesh::countFacesWithSides(size_t n) const {
    size_t count = 0;
    for (const auto& f : faces)
        if (f.size() == n) ++count;
    return count;
}

Mesh readObj(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open " + path);

    Mesh mesh;
    std::string line;
    size_t lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (line.size() < 2) continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "v") {
            std::array<double, 3> p{};
            if (!(ss >> p[0] >> p[1] >> p[2]))
                throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": bad vertex");
            mesh.vertices.push_back(p);
        } else if (tag == "f") {
            std::vector<int> face;
            std::string token;
            while (ss >> token) {
                // Accept v, v/t, v/t/n and v//n; only the position index matters here.
                long idx = std::strtol(token.c_str(), nullptr, 10);
                if (idx == 0)
                    throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": bad face index");
                long resolved = idx > 0 ? idx - 1 : static_cast<long>(mesh.vertices.size()) + idx;
                if (resolved < 0 || resolved >= static_cast<long>(mesh.vertices.size()))
                    throw std::runtime_error(path + ":" + std::to_string(lineNo) + ": face index out of range");
                face.push_back(static_cast<int>(resolved));
            }
            if (face.size() >= 3) mesh.faces.push_back(std::move(face));
        }
    }
    if (mesh.vertices.empty() || mesh.faces.empty())
        throw std::runtime_error(path + ": mesh has no geometry");
    return mesh;
}

void writeObj(const std::string& path, const Mesh& mesh) {
    FILE* out = std::fopen(path.c_str(), "w");
    if (!out) throw std::runtime_error("cannot write " + path);
    std::fprintf(out, "# remersher\n");
    for (const auto& v : mesh.vertices) std::fprintf(out, "v %.9g %.9g %.9g\n", v[0], v[1], v[2]);
    for (const auto& f : mesh.faces) {
        std::fputc('f', out);
        for (int idx : f) std::fprintf(out, " %d", idx + 1);
        std::fputc('\n', out);
    }
    if (std::fclose(out) != 0) throw std::runtime_error("error writing " + path);
}

}  // namespace remersher
