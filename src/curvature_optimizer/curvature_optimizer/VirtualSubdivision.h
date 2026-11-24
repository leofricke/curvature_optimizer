#pragma once

#include "CurvatureOptimizerWorkspace.h"
#include "utils/Barycentric.h"
#include "utils/InlineVector.h"
#include <pmp/surface_mesh.h>

struct VirtualTriangle;

/**
* Subdivides a mesh face into virtual triangles.
*/
struct VirtualSubdivision
{
    VirtualSubdivision() = default;

    void updateCenterWeights(const CurvatureOptimizerWorkspace& ws, pmp::Face f);
    void updateCenterWeights(const pmp::SurfaceMesh& mesh, pmp::Face f);
    Eigen::Vector3d getCenterPoint(const CurvatureOptimizerWorkspace& ws, pmp::Face f) const;
    Eigen::Vector3d getCenterPoint(const pmp::SurfaceMesh& mesh, pmp::Face f) const;
    double getVertexWeight(const CurvatureOptimizerWorkspace& ws, pmp::Face f, pmp::Vertex v) const;

    void initVirtualTriangles(const pmp::SurfaceMesh& mesh, pmp::Face f, InlineVector<VirtualTriangle, 4>& triangles) const;

    InlineVector<double, 4> m_centerWeights;

private:
    void updateCenterWeights(const Eigen::Matrix<double, Eigen::Dynamic, 3>& poly);
};

/**
* A virtual triangle that is subdivided into six more virtual triangles by splitting at the incircle center.
* Represents the initial (target) triangle shape that is preserved by the conformal energy.
*/
struct VirtualTriangle
{
    VirtualTriangle(pmp::Halfedge h, const Eigen::Vector3d& p0, const Eigen::Vector3d& p1, const Eigen::Vector3d& center);

    pmp::Halfedge m_halfedge;
    double m_le01;
    double m_le12;
    double m_le20;
    double m_area;
    BarycentricCoords m_incenter;
    BarycentricCoords m_q01;
    BarycentricCoords m_q12;
    BarycentricCoords m_q20;
};
