#pragma once

#include "VirtualSubdivision.h"
#include "utils/Barycentric.h"
#include "utils/InlineVector.h"

#include <Eigen/Core>
#include <pmp/surface_mesh.h>

struct TargetVertexCurvature
{
    pmp::Vertex m_vertex;
    double m_kmax;
    double m_kmin;
};

struct TargetBoundaryCurvature
{
    pmp::Vertex m_vertex;
    double m_k;
};

struct TargetPosition
{
    BarycentricPoint m_bary;
    Eigen::Vector3d m_position;
};

struct TargetMetric
{
    pmp::Face m_face;
    double m_area;
    VirtualSubdivision m_subdivision;
    InlineVector<VirtualTriangle, 4> m_triangles;
};

TargetMetric createTargetMetric(const pmp::SurfaceMesh& mesh, pmp::Face f);
