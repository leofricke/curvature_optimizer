#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

double computeIntegratedBoundaryCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v);
double computeBoundaryCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v);

double computeIntegratedBoundaryCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d>& g);
double computeBoundaryCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d>& g);
