#include "CommonDerivatives.h"

#include <Eigen/Dense>
#include "utils/EigenMathUtil.h"

void computeCrossProductJacobian(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, Eigen::Matrix<double, 3, 9>& J)
{
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e20 = p0 - p2;
    J.block<3, 3>(0, 3) = crossMatrix<double>(e20);
    J.block<3, 3>(0, 6) = crossMatrix<double>(e01);
    J.block<3, 3>(0, 0) = -J.block<3, 3>(0, 3) - J.block<3, 3>(0, 6);
}

void computeTriangleNormalJacobian(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, Eigen::Matrix<double, 3, 9>& J)
{
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e20 = p0 - p2;
    Eigen::Vector3d n = e01.cross(-e20);
    const double doubleArea = n.norm();
    n.normalize();
    const Eigen::Matrix3d P = (Eigen::Matrix3d::Identity() - n * n.transpose()) / doubleArea;

    J.block<3, 3>(0, 3) = P * crossMatrix<double>(e20);
    J.block<3, 3>(0, 6) = P * crossMatrix<double>(e01);
    J.block<3, 3>(0, 0) = -J.block<3, 3>(0, 3) - J.block<3, 3>(0, 6);
}

void computeCotanGradient(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, Eigen::RowVector<double, 9>& g)
{
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e20 = p0 - p2;
    Eigen::Vector3d n = e01.cross(-e20);
    const double doubleArea = n.norm();
    n.normalize();
    const double d = e01.dot(-e20);

    // triangle area derivs (times two)
    const Eigen::Vector3d dav1 = n.cross(e20);
    const Eigen::Vector3d dav2 = n.cross(e01);

    g.segment<3>(3) = (-e20 + d * -dav1 / doubleArea) / doubleArea;
    g.segment<3>(6) = (e01 + d * -dav2 / doubleArea) / doubleArea;
    g.segment<3>(0) = -g.segment<3>(3) - g.segment<3>(6);
}

void computeNormalizedOuterProductJacobian(const Eigen::Vector3d& e, std::array<Eigen::Matrix3d, 3>& J)
{
    const double le = e.norm();
    const Eigen::Vector3d e0 = e.normalized();

    const Eigen::Matrix3d ee0 = e0 * e0.transpose();
    const std::array<Eigen::Vector3d, 3> units = { Eigen::Vector3d { 1.0, 0.0, 0.0 }, Eigen::Vector3d { 0.0, 1.0, 0.0 }, Eigen::Vector3d { 0.0, 0.0, 1.0 } };

    for (size_t i = 0; i < 3; ++i)
    {
        J[i] = (2.0 * sym(e0 * units[i].transpose()) - ee0 * 2.0 * e0[i]) / le;
    }
}
