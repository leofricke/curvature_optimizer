#pragma once

#include "IVertexCurvatureOperator.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

struct MeanGaussianVertexCurvature : public IVertexCurvatureOperator
{
    MeanGaussianVertexCurvature();

    Eigen::Vector2d computePrincipalCurvatures(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v) override;
    Eigen::Vector2d computePrincipalCurvaturesGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& gkmax,
                                                       LocalGradient<Eigen::RowVector3d, 8>& gkmin) override;
};
