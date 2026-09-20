#include "WritePLY.h"

#include "happly.h"

void writeMeshPLY(const std::filesystem::path& file, const pmp::SurfaceMesh& mesh, const std::unordered_map<std::string, std::vector<double>>& vertexData)
{
    std::vector<std::array<double, 3>> meshVertexPositions;
    std::vector<std::vector<size_t>> meshFaceIndices;

    meshVertexPositions.reserve(mesh.n_vertices());
    meshFaceIndices.reserve(mesh.n_faces());
    for (pmp::Vertex v : mesh.vertices())
    {
        const auto& pos = mesh.position(v);
        meshVertexPositions.push_back({ pos[0], pos[1], pos[2] });
    }
    for (pmp::Face f : mesh.faces())
    {
        std::vector<size_t> faceIndices;
        for (pmp::Vertex v : mesh.vertices(f))
        {
            faceIndices.push_back(v.idx());
        }
        meshFaceIndices.emplace_back(std::move(faceIndices));
    }

    happly::PLYData ply;
    ply.addVertexPositions(meshVertexPositions);
    ply.addFaceIndices(meshFaceIndices);

    for (const auto& [name, data] : vertexData)
    {
        ply.getElement("vertex").addProperty<double>(name, data);
    }

    ply.write(file.string(), happly::DataFormat::ASCII);
}
