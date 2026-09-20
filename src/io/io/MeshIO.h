#pragma once

#include <pmp/surface_mesh.h>
#include <filesystem>

void readMesh(std::filesystem::path file, pmp::SurfaceMesh& mesh);
void writeMesh(std::filesystem::path file, const pmp::SurfaceMesh& mesh);
