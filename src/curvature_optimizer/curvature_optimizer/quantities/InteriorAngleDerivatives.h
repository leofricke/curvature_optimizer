#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

double computeInteriorAngle(const CurvatureOptimizerWorkspace& ws, pmp::Halfedge h);
double computeInteriorAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Halfedge h, LocalGradient<Eigen::RowVector3d, 3>& g);
