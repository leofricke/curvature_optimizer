#pragma once

#include <filesystem>
#include <unordered_map>
#include <vector>
#include <string>
#include <pmp/surface_mesh.h>

void writeMeshPLY(const std::filesystem::path& file, const pmp::SurfaceMesh& mesh, const std::unordered_map<std::string, std::vector<double>>& vertexData = {});
