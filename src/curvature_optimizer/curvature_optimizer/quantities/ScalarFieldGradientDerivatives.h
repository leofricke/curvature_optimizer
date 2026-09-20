#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

struct ScalarFieldFaceData
{
    pmp::Face m_face;
    double m_f0;
    double m_f1;
    double m_f2;
    const LocalGradient<Eigen::RowVector3d, 8>& m_g0;
    const LocalGradient<Eigen::RowVector3d, 8>& m_g1;
    const LocalGradient<Eigen::RowVector3d, 8>& m_g2;
};

Eigen::Vector3d computeScalarFieldGradient(const CurvatureOptimizerWorkspace& ws, const ScalarFieldFaceData& fieldData);
Eigen::Vector3d computeScalarFieldGradientJacobian(const CurvatureOptimizerWorkspace& ws, const ScalarFieldFaceData& fieldData, LocalGradient<Eigen::Matrix3d, 20>& J);
