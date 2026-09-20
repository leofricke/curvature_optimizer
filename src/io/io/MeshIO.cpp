#include "MeshIO.h"

#include "ReadPLY.h"
#include "WritePLY.h"

#include <pmp/io/io.h>

void readMesh(std::filesystem::path file, pmp::SurfaceMesh& mesh)
{
    mesh.clear();

    file.make_preferred();
    const std::string ext = file.extension().string();
    if (ext == ".ply")
    {
        readMeshPLY(file, mesh);
    }
    else if (ext == ".obj" || ext == ".off" || ext == ".pmp" || ext == ".stl")
    {
        read(mesh, file);
    }
    else
    {
        throw std::runtime_error("Unsupported file extension");
    }
}

void writeMesh(std::filesystem::path file, const pmp::SurfaceMesh& mesh)
{
    file.make_preferred();
    const std::string ext = file.extension().string();
    if (ext == ".ply")
    {
        writeMeshPLY(file, mesh);
    }
    else if (ext == ".obj" || ext == ".off" || ext == ".pmp" || ext == ".stl")
    {
        write(mesh, file);
    }
    else
    {
        throw std::runtime_error("Unsupported file extension");
    }
}