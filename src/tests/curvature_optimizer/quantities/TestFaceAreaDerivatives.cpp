#include <gtest/gtest.h>

#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/quantities/FaceAreaDerivatives.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class TriangleAreaDerivativesTest : public testing::Test
{
public:
    TriangleAreaDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(TriangleAreaDerivativesTest, ComputeTriangleAreaGradient)
{
    const double eps = 1e-8;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Face f : mesh.faces())
    {
        LocalGradient<Eigen::RowVector3d> g;
        computeTriangleAreaGradient(ws, f, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double A_plus = computeTriangleArea(ws, f);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double A_minus = computeTriangleArea(ws, f);

                const double dA = (A_plus - A_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], dA));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}

TEST_F(TriangleAreaDerivativesTest, ComputeTriangleNormalJacobian)
{
    const double eps = 1e-7;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Face f : mesh.faces())
    {
        LocalGradient<Eigen::Matrix3d> g;
        computeTriangleNormalJacobian(ws, f, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const Eigen::Vector3d n_plus = computeTriangleNormal(ws, f);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const Eigen::Vector3d n_minus = computeTriangleNormal(ws, f);

                const Eigen::Vector3d dn = (n_plus - n_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDerivMatrix(g[vv].col(j), dn));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}

class PolygonAreaDerivativesTest : public testing::Test
{
public:
    PolygonAreaDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("two_quads.obj"));

        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(PolygonAreaDerivativesTest, ComputePolygonAreaGradient)
{
    const double eps = 1e-8;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Face f : mesh.faces())
    {
        LocalGradient<Eigen::RowVector3d> g;
        computePolygonAreaGradient(ws, f, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double A_plus = computePolygonArea(ws, f);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double A_minus = computePolygonArea(ws, f);

                const double dA = (A_plus - A_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], dA));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}

TEST_F(PolygonAreaDerivativesTest, ComputePolygonNormalJacobian)
{
    const double eps = 1e-7;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Face f : mesh.faces())
    {
        LocalGradient<Eigen::Matrix3d> g;
        computePolygonNormalJacobian(ws, f, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const Eigen::Vector3d n_plus = computePolygonNormal(ws, f);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const Eigen::Vector3d n_minus = computePolygonNormal(ws, f);

                const Eigen::Vector3d dn = (n_plus - n_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDerivMatrix(g[vv].col(j), dn));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}