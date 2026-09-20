#include <gtest/gtest.h>

#include "curvature_optimizer/quantities/CommonDerivatives.h"
#include "DerivativeTestUtil.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class CommonDerivativesTest : public ::testing::Test
{
public:
    CommonDerivativesTest()
    {
        read(mesh, getTestFile("vertex_onering.off"));
    }

    pmp::SurfaceMesh mesh;
};

TEST_F(CommonDerivativesTest, ComputeCrossProductJacobian)
{
    const auto func = [&](const std::array<Eigen::Vector3d, 3>& args)
    {
        return (args[1] - args[0]).cross(args[2] - args[0]);
    };

    for (pmp::Face f : mesh.faces())
    {
        const pmp::Halfedge h = mesh.halfedge(f);
        const Eigen::Vector3d p0 = mesh.position(mesh.from_vertex(h));
        const Eigen::Vector3d p1 = mesh.position(mesh.to_vertex(h));
        const Eigen::Vector3d p2 = mesh.position(mesh.to_vertex(mesh.next_halfedge(h)));

        Eigen::Matrix<double, 3, 9> J;
        computeCrossProductJacobian(p0, p1, p2, J);

        std::array<std::array<Eigen::Vector3d, 3>, 3> derivs;
        numericalDeriv<Eigen::Vector3d, Eigen::Vector3d, 3>(func, { p0, p1, p2 }, derivs);

        for (size_t i = 0; i < 9; ++i)
        {
            EXPECT_TRUE(checkDerivMatrix(J.col(i), derivs[i / 3][i % 3]));
        }
    }
}

TEST_F(CommonDerivativesTest, ComputeTriangleNormalJacobian)
{
    const auto func = [&](const std::array<Eigen::Vector3d, 3>& args)
    {
        return (args[1] - args[0]).cross(args[2] - args[0]).normalized();
    };

    for (pmp::Face f : mesh.faces())
    {
        const pmp::Halfedge h = mesh.halfedge(f);
        const Eigen::Vector3d p0 = mesh.position(mesh.from_vertex(h));
        const Eigen::Vector3d p1 = mesh.position(mesh.to_vertex(h));
        const Eigen::Vector3d p2 = mesh.position(mesh.to_vertex(mesh.next_halfedge(h)));

        Eigen::Matrix<double, 3, 9> J;
        computeTriangleNormalJacobian(p0, p1, p2, J);

        std::array<std::array<Eigen::Vector3d, 3>, 3> derivs;
        numericalDeriv<Eigen::Vector3d, Eigen::Vector3d, 3>(func, { p0, p1, p2 }, derivs);

        for (size_t i = 0; i < 9; ++i)
        {
            EXPECT_TRUE(checkDerivMatrix(J.col(i), derivs[i / 3][i % 3]));
        }
    }
}

TEST_F(CommonDerivativesTest, ComputeCotanGradient)
{
    const auto func = [&](const std::array<Eigen::Vector3d, 3>& args)
    {
        return (args[1] - args[0]).dot(args[2] - args[0]) / (args[1] - args[0]).cross(args[2] - args[0]).norm();
    };

    for (pmp::Face f : mesh.faces())
    {
        const pmp::Halfedge h = mesh.halfedge(f);
        const Eigen::Vector3d p0 = mesh.position(mesh.from_vertex(h));
        const Eigen::Vector3d p1 = mesh.position(mesh.to_vertex(h));
        const Eigen::Vector3d p2 = mesh.position(mesh.to_vertex(mesh.next_halfedge(h)));

        Eigen::RowVector<double, 9> g;
        computeCotanGradient(p0, p1, p2, g);

        std::array<std::array<double, 3>, 3> derivs;
        numericalDeriv<Eigen::Vector3d, double, 3>(func, { p0, p1, p2 }, derivs);

        for (size_t i = 0; i < 9; ++i)
        {
            EXPECT_TRUE(checkDeriv(g[i], derivs[i / 3][i % 3]));
        }
    }
}

TEST(CommonDerivativesTest2, ComputeNormalizedOuterProductJacobian)
{
    const auto func = [&](const std::array<Eigen::Vector3d, 2>& args)
    {
        const Eigen::Vector3d e0 = (args[1] - args[0]).normalized();
        const Eigen::Matrix3d ee = e0 * e0.transpose();
        return ee;
    };

    const Eigen::Vector3d p0 = { 0.0, 0.0, 0.0 };
    const Eigen::Vector3d p1 = { 1.0, 2.0, 3.0 };
    const Eigen::Vector3d e = p1 - p0;

    std::array<Eigen::Matrix3d, 3> J;
    computeNormalizedOuterProductJacobian(e, J);

    std::array<std::array<Eigen::Matrix3d, 3>, 2> derivs;
    numericalDeriv<Eigen::Vector3d, Eigen::Matrix3d, 2>(func, { p0, p1 }, derivs);

    EXPECT_TRUE(checkDerivMatrix(-J[0], derivs[0][0]));
    EXPECT_TRUE(checkDerivMatrix(-J[1], derivs[0][1]));
    EXPECT_TRUE(checkDerivMatrix(-J[2], derivs[0][2]));
    EXPECT_TRUE(checkDerivMatrix(J[0], derivs[1][0]));
    EXPECT_TRUE(checkDerivMatrix(J[1], derivs[1][1]));
    EXPECT_TRUE(checkDerivMatrix(J[2], derivs[1][2]));
}