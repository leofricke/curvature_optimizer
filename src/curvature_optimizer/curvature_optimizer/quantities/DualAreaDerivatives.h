#pragma once

#include "LocalGradient.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

enum class DualAreaMode
{
    Barycentric = 0,
    Voronoi = 1,
    VoronoiMixed = 2
};
static std::array<std::string, 3> s_dualAreaModeNames = { "Barycentric", "Voronoi", "VoronoiMixed" };

double computeDualArea(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, DualAreaMode mode = DualAreaMode::Barycentric);
double computeDualAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g, DualAreaMode mode = DualAreaMode::Barycentric);

double computeBarycentricArea(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v);
double computeBarycentricAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g);

double computeVoronoiArea(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, bool mixed = false);
double computeVoronoiAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g, bool mixed = false);
