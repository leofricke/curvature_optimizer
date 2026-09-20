#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

double computeVertexIntegratedMeanCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v);
double computeVertexMeanCurvature(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v);

double computeVertexIntegratedMeanCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g);
double computeVertexMeanCurvatureGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g);
