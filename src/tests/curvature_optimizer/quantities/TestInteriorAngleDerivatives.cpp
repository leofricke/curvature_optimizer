#include <gtest/gtest.h>

#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "curvature_optimizer/quantities/InteriorAngleDerivatives.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class InteriorAngleDerivativesTest : public testing::Test
{
public:
    InteriorAngleDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(InteriorAngleDerivativesTest, ComputeInteriorAngleGradient)
{
    const double eps = 1e-8;
    const double atol = 1e-7;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Halfedge h : mesh.halfedges())
    {
        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::RowVector3d, 3> g;
        computeInteriorAngleGradient(ws, h, g);

        for (pmp::Vertex vv : g)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const double a_plus = computeInteriorAngle(ws, h);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const double a_minus = computeInteriorAngle(ws, h);

                const double da = (a_plus - a_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDeriv(g[vv][j], da, atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}
