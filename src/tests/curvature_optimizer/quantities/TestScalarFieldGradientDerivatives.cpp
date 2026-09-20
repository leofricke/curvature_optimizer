#include <gtest/gtest.h>

#include "curvature_optimizer/quantities/ScalarFieldGradientDerivatives.h"
#include "curvature_optimizer/CurvatureOptimizerWorkspace.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class ScalarFieldGradientDerivativesTest : public testing::Test
{
public:
    ScalarFieldGradientDerivativesTest()
        : ws(mesh)
    {
        read(mesh, getTestFile("cylinder_low.off"));

        ws.initMeshData();
        ws.m_dofs.initMeshDofs({});
    }

    pmp::SurfaceMesh mesh;
    CurvatureOptimizerWorkspace ws;
};

TEST_F(ScalarFieldGradientDerivativesTest, ComputeScalarFieldGradientJacobian)
{
    const double eps = 1e-5;
    const double atol = 1e-7;
    Eigen::VectorXd& dofs = ws.m_dofs.getDofs();
    for (pmp::Face f : mesh.faces())
    {
        const ScalarFieldFaceData fieldData = { .m_face = f, .m_f0 = 2.0, .m_f1 = 3.0, .m_f2 = 4.0, .m_g0 = {}, .m_g1 = {}, .m_g2 = {} };

        ws.updateGeometry(EnergyFunction::DerivativeMode::FirstDerivs);
        LocalGradient<Eigen::Matrix3d, 20> J;
        computeScalarFieldGradientJacobian(ws, fieldData, J);

        for (pmp::Vertex vv : J)
        {
            for (size_t j = 0; j < 3; ++j)
            {
                const size_t i = ws.m_dofs.dofIdx(vv) * 3 + j;
                const double orig = dofs[i];

                dofs[i] = orig + eps;
                ws.updateGeometry();
                const Eigen::Vector3d g_plus = computeScalarFieldGradient(ws, fieldData);

                dofs[i] = orig - eps;
                ws.updateGeometry();
                const Eigen::Vector3d g_minus = computeScalarFieldGradient(ws, fieldData);

                const Eigen::Vector3d dg = (g_plus - g_minus) / (2.0 * eps);
                EXPECT_TRUE(checkDerivMatrix(J[vv].col(j), dg, atol));
                dofs[i] = orig;
                ws.updateGeometry();
            }
        }
    }
}
