#pragma once

#include <Eigen/Core>

void eigendecomposition3(const Eigen::Matrix3d& A, Eigen::Vector3d& w);
void eigendecomposition3(const Eigen::Matrix3d& A, Eigen::Vector3d& w, Eigen::Matrix3d& Q);
