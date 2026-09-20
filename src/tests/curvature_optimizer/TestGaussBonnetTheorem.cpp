#include <gtest/gtest.h>

#include "curvature_optimizer/quantities/BoundaryCurvatureDerivatives.h"
#include "curvature_optimizer/quantities/GaussianCurvatureDerivatives.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>
#include <numbers>

TEST(GaussBonnetTheorem, GaussBonnetVertexSphere)
{
    pmp::SurfaceMesh mesh;
    read(mesh, getTestFile("sphere.off"));
    CurvatureOptimizerWorkspace ws(mesh);
    ws.initMeshData();

    double sum = 0.0;
    for (pmp::Vertex v : mesh.vertices())
    {
        sum += computeVertexIntegratedGaussianCurvature(ws, v);
    }

    EXPECT_NEAR(sum, 4.0 * std::numbers::pi, 1e-7);
}

TEST(GaussBonnetTheorem, GaussBonnetVertexHemisphere)
{
    pmp::SurfaceMesh mesh;
    read(mesh, getTestFile("hemisphere.off"));
    CurvatureOptimizerWorkspace ws(mesh);
    ws.initMeshData();

    double sum = 0.0;
    for (pmp::Vertex v : mesh.vertices())
    {
        if (!mesh.is_boundary(v))
        {
            sum += computeVertexIntegratedGaussianCurvature(ws, v);
        }
        else
        {
            sum += computeIntegratedBoundaryCurvature(ws, v);
        }
    }

    EXPECT_NEAR(sum, 2.0 * std::numbers::pi, 1e-7);
}
