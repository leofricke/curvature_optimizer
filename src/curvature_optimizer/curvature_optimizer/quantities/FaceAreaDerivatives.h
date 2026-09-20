#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

double computeFaceArea(const CurvatureOptimizerWorkspace& ws, pmp::Face f);
double computeFaceAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::RowVector3d, 4>& g);

double computeTriangleArea(const CurvatureOptimizerWorkspace& ws, pmp::Face f);
double computeTriangleAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::RowVector3d, 4>& g);

double computePolygonArea(const CurvatureOptimizerWorkspace& ws, pmp::Face f);
double computePolygonAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::RowVector3d, 4>& g);

Eigen::Vector3d computeFaceNormal(const CurvatureOptimizerWorkspace& ws, pmp::Face f);
Eigen::Vector3d computeFaceNormalJacobian(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::Matrix3d, 4>& J);

Eigen::Vector3d computeTriangleNormal(const CurvatureOptimizerWorkspace& ws, pmp::Face f);
Eigen::Vector3d computeTriangleNormalJacobian(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::Matrix3d, 4>& J);

Eigen::Vector3d computePolygonNormal(const CurvatureOptimizerWorkspace& ws, pmp::Face f);
Eigen::Vector3d computePolygonNormalJacobian(const CurvatureOptimizerWorkspace& ws, pmp::Face f, LocalGradient<Eigen::Matrix3d, 4>& J);
