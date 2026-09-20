#pragma once

#include <filesystem>
#include <pmp/surface_mesh.h>

void readMeshPLY(const std::filesystem::path& file, pmp::SurfaceMesh& mesh);
