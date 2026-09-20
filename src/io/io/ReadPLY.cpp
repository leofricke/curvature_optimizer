#include "ReadPLY.h"

#include "happly.h"

void readMeshPLY(const std::filesystem::path& file, pmp::SurfaceMesh& mesh)
{
    happly::PLYData ply(file.string());

    const std::vector<std::array<double, 3>> meshVertexPositions = ply.getVertexPositions();
    const std::vector<std::vector<size_t>> meshFaceIndices = ply.getFaceIndices<size_t>();

    mesh.clear();
    for (const auto& pos : meshVertexPositions)
    {
        mesh.add_vertex(pmp::Point(pos[0], pos[1], pos[2]));
    }
    for (const auto& face : meshFaceIndices)
    {
        std::vector<pmp::Vertex> vertices;
        for (const auto& vertexIndex : face)
        {
            vertices.push_back(pmp::Vertex(vertexIndex));
        }
        mesh.add_face(vertices);
    }
}
