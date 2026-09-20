#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

double computeVertexIntegratedGaussianCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v);
double computeVertexGaussianCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v);

double computeVertexIntegratedGaussianCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g);
double computeVertexGaussianCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g);
