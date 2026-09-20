#include <gtest/gtest.h>

#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/quantities/DihedralAngleDerivatives.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class TriangleDihedralAngleDerivativesTest : public testing::Test
{
public:
    TriangleDihedralAngleDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(TriangleDihedralAngleDerivativesTest, ComputeTriangleDihedralAngleGradient)
{
    const double eps = 1e-8;
    const double atol = 1e-6;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Edge e : mesh.edges())
    {
        LocalGradient<Eigen::RowVector3d> g;
        computeTriangleDihedralAngleGradient(ws, e, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double b_plus = computeTriangleDihedralAngle(ws, e);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double b_minus = computeTriangleDihedralAngle(ws, e);

                const double db = (b_plus - b_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], db, atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}

class PolygonDihedralAngleDerivativesTest : public testing::Test
{
public:
    PolygonDihedralAngleDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("two_quads_bend.obj"));

        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(PolygonDihedralAngleDerivativesTest, ComputePolygonDihedralAngleGradient)
{
    const double eps = 1e-8;
    const double atol = 1e-6;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Edge e : mesh.edges())
    {
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::RowVector3d> g;
        computePolygonDihedralAngleGradient(ws, e, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double b_plus = computePolygonDihedralAngle(ws, e);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double b_minus = computePolygonDihedralAngle(ws, e);

                const double db = (b_plus - b_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], db, atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}
