#include "MeshSequenceExporter.h"

#include "WritePLY.h"

MeshSequenceExporter::MeshSequenceExporter(const std::filesystem::path& outputDir)
    : m_outputDir(outputDir), m_counter(0)
{
}

void MeshSequenceExporter::reset()
{
    m_counter = 0;

    if (std::filesystem::exists(m_outputDir))
    {
        for (const auto& entry : std::filesystem::directory_iterator(m_outputDir))
        {
            if (entry.is_regular_file())
            {
                const std::string name = entry.path().filename().string();
                if (name.size() > 4 && name.substr(0, 5) == "mesh_" && name.substr(name.size() - 4) == ".ply")
                {
                    std::filesystem::remove(entry.path());
                }
            }
        }
    }
}

void MeshSequenceExporter::appendMesh(const pmp::SurfaceMesh& mesh, const std::unordered_map<std::string, std::vector<double>>& vertexData)
{
    if (!std::filesystem::exists(m_outputDir))
    {
        std::filesystem::create_directories(m_outputDir);
    }

    std::string s = std::to_string(m_counter);
    const size_t width = 4;
    if (s.size() < width)
    {
        s.insert(0, width - s.size(), '0');
    }
    std::filesystem::path outputFile = m_outputDir / ("mesh_" + s + ".ply");

    writeMeshPLY(outputFile, mesh, vertexData);
    m_counter++;
}
