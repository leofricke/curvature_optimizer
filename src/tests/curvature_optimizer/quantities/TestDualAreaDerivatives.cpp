#include <gtest/gtest.h>

#include "curvature_optimizer/quantities/DualAreaDerivatives.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class DualAreaDerivativesTest : public testing::Test
{
public:
    DualAreaDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("vertex_onering.off"));

        ws.m_dualAreaMode = DualAreaMode::Barycentric;
        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(DualAreaDerivativesTest, ComputeBarycentricAreaGradient)
{
    const double eps = 1e-8;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Vertex v0 : mesh.vertices())
    {
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::RowVector3d, 8> g;
        computeBarycentricAreaGradient(ws, v0, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double A_plus = computeBarycentricArea(ws, v0);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double A_minus = computeBarycentricArea(ws, v0);

                const double dA = (A_plus - A_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], dA));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}

TEST_F(DualAreaDerivativesTest, ComputeVoronoiAreaGradient)
{
    const double eps = 1e-8;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Vertex v0 : mesh.vertices())
    {
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::RowVector3d, 8> g;
        computeVoronoiAreaGradient(ws, v0, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double A_plus = computeVoronoiArea(ws, v0);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double A_minus = computeVoronoiArea(ws, v0);

                const double dA = (A_plus - A_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], dA));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}

TEST(VoronoiArea, EquilateralTriangle)
{
    pmp::SurfaceMesh mesh;
    const pmp::Vertex v0 = mesh.add_vertex(pmp::Point(0.0, 0.0, 0.0));
    const pmp::Vertex v1 = mesh.add_vertex(pmp::Point(1.0, 0.0, 0.0));
    const pmp::Vertex v2 = mesh.add_vertex(pmp::Point(0.5, std::sqrt(3.0) * 0.5, 0.0));
    mesh.add_triangle(v0, v1, v2);
    CurvatureOptimizerWorkspace ws(mesh);
    ws.initMeshData();

    const double a0 = computeVoronoiArea(ws, v0);
    const double a1 = computeVoronoiArea(ws, v1);
    const double a2 = computeVoronoiArea(ws, v2);
    const double expected = (std::sqrt(3.0) / 4.0) / 3.0;
    EXPECT_NEAR(a0, expected, 1e-8);
    EXPECT_NEAR(a1, expected, 1e-8);
    EXPECT_NEAR(a2, expected, 1e-8);
}

TEST(VoronoiArea, RightTriangle)
{
    pmp::SurfaceMesh mesh;
    const pmp::Vertex v0 = mesh.add_vertex(pmp::Point(0.0, 0.0, 0.0));
    const pmp::Vertex v1 = mesh.add_vertex(pmp::Point(1.0, 0.0, 0.0));
    const pmp::Vertex v2 = mesh.add_vertex(pmp::Point(0.0, 1.0, 0.0));
    mesh.add_triangle(v0, v1, v2);
    CurvatureOptimizerWorkspace ws(mesh);
    ws.initMeshData();

    const double a0 = computeVoronoiArea(ws, v0);
    const double a1 = computeVoronoiArea(ws, v1);
    const double a2 = computeVoronoiArea(ws, v2);
    EXPECT_NEAR(a0, 0.25, 1e-12);
    EXPECT_NEAR(a1, 0.125, 1e-12);
    EXPECT_NEAR(a2, 0.125, 1e-12);
}

TEST(VoronoiArea, ObtuseTriangle)
{
    pmp::SurfaceMesh mesh;
    const pmp::Vertex v0 = mesh.add_vertex(pmp::Point(0.25, 0.25, 0.0));
    const pmp::Vertex v1 = mesh.add_vertex(pmp::Point(1.0, 0.0, 0.0));
    const pmp::Vertex v2 = mesh.add_vertex(pmp::Point(0.0, 1.0, 0.0));
    mesh.add_triangle(v0, v1, v2);
    CurvatureOptimizerWorkspace ws(mesh);
    ws.initMeshData();

    const double a0 = computeVoronoiArea(ws, v0);
    const double a1 = computeVoronoiArea(ws, v1);
    const double a2 = computeVoronoiArea(ws, v2);
    EXPECT_NEAR(a0 + a1 + a2, 0.25, 1e-12);
    EXPECT_GT(a0, 0.0);
    EXPECT_LT(a1, 0.0);
    EXPECT_LT(a2, 0.0);
}
