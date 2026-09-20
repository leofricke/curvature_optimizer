#pragma once

#include <filesystem>
#include <pmp/surface_mesh.h>

class MeshSequenceExporter
{
public:
    MeshSequenceExporter(const std::filesystem::path& outputDir);

    void reset();
    void appendMesh(const pmp::SurfaceMesh& mesh, const std::unordered_map<std::string, std::vector<double>>& vertexData = {});

private:
    std::filesystem::path m_outputDir;
    size_t m_counter;
};
