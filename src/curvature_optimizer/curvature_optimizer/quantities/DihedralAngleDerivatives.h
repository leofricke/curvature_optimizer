#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

double computeDihedralAngle(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, bool modifyAngle = false);
double computeDihedralAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, LocalGradient<Eigen::RowVector3d, 4>& g, bool modifyAngle = false);

double computeTriangleDihedralAngle(const CurvatureOptimizerWorkspace& ws, pmp::Edge e);
double computeTriangleDihedralAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, LocalGradient<Eigen::RowVector3d, 4>& g);

double computePolygonDihedralAngle(const CurvatureOptimizerWorkspace& ws, pmp::Edge e);
double computePolygonDihedralAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, LocalGradient<Eigen::RowVector3d, 4>& g);
