#pragma once

#include <array>
#include <Eigen/Core>

void computeCrossProductJacobian(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, Eigen::Matrix<double, 3, 9>& J);
void computeTriangleNormalJacobian(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, Eigen::Matrix<double, 3, 9>& J);
void computeCotanGradient(const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, Eigen::RowVector<double, 9>& g);
void computeNormalizedOuterProductJacobian(const Eigen::Vector3d& e, std::array<Eigen::Matrix3d, 3>& J);
