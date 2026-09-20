#pragma once

#include "PrincipalCurvatures.h"
#include "quantities/LocalGradient.h"

#include <array>
#include <Eigen/Core>
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

enum class VertexCurvatureOperator
{
    CurvatureTensor = 0,
    MeanGaussian = 1
};
static std::array<std::string, 2> s_vertexCurvatureOperatorNames = { "CurvatureTensor", "MeanGaussian" };

struct IVertexCurvatureOperator
{
    virtual ~IVertexCurvatureOperator() = default;

    virtual Eigen::Vector2d computePrincipalCurvatures(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v) = 0;
    virtual Eigen::Vector2d computePrincipalCurvaturesGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& gkmax,
                                                               LocalGradient<Eigen::RowVector3d, 8>& gkmin) = 0;

    PrincipalCurvatures m_curv;
};
