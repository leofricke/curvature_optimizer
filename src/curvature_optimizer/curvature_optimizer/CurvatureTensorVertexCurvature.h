#pragma once

#include "IVertexCurvatureOperator.h"

#include <Eigen/Core>
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

struct CurvatureTensorVertexCurvature : public IVertexCurvatureOperator
{
    CurvatureTensorVertexCurvature();

    Eigen::Matrix3d computeTensor(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v) const;
    Eigen::Matrix3d computeTensorDeriv(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::Matrix<double, 3, 9>>& dT) const;

    Eigen::Vector2d computePrincipalCurvatures(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v) override;
    Eigen::Vector2d computePrincipalCurvaturesGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& gkmax,
                                                       LocalGradient<Eigen::RowVector3d, 8>& gkmin) override;

    void updateCurvaturesAndDirections(const Eigen::Matrix3d& tensor);
};
