#include "ScalarFieldGradientDerivatives.h"

#include "DiffRules.h"
#include "../CurvatureOptimizerWorkspace.h"
#include "utils/EigenMathUtil.h"

Eigen::Vector3d computeScalarFieldGradient(const CurvatureOptimizerWorkspace& ws, const ScalarFieldFaceData& fieldData)
{
    if (ws.m_mesh.valence(fieldData.m_face) != 3)
    {
        throw std::runtime_error("Scalar Field Gradient only supported on triangles");
    }

    Eigen::Vector3d grad = Eigen::Vector3d::Zero();

    const pmp::Halfedge h = ws.m_mesh.halfedge(fieldData.m_face);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e12 = p2 - p1;
    const Eigen::Vector3d e20 = p0 - p2;

    Eigen::Vector3d n = e01.cross(-e20);
    const double doubleArea = n.norm();
    n.normalize();
    const Eigen::Vector3d g0 = n.cross(e12) / doubleArea;
    const Eigen::Vector3d g1 = n.cross(e20) / doubleArea;
    const Eigen::Vector3d g2 = n.cross(e01) / doubleArea;
    grad += fieldData.m_f0 * g0;
    grad += fieldData.m_f1 * g1;
    grad += fieldData.m_f2 * g2;

    return grad;
}

Eigen::Vector3d computeScalarFieldGradientJacobian(const CurvatureOptimizerWorkspace& ws, const ScalarFieldFaceData& fieldData, LocalGradient<Eigen::Matrix3d, 20>& J)
{
    if (ws.m_mesh.valence(fieldData.m_face) != 3)
    {
        throw std::runtime_error("Scalar Field Gradient only supported on triangles");
    }

    Eigen::Vector3d grad = Eigen::Vector3d::Zero();

    const pmp::Halfedge h = ws.m_mesh.halfedge(fieldData.m_face);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e12 = p2 - p1;
    const Eigen::Vector3d e20 = p0 - p2;

    Eigen::Vector3d n = e01.cross(-e20);
    const double doubleArea = n.norm();
    n.normalize();
    const Eigen::Vector3d g0 = n.cross(e12) / doubleArea;
    const Eigen::Vector3d g1 = n.cross(e20) / doubleArea;
    const Eigen::Vector3d g2 = n.cross(e01) / doubleArea;
    grad += fieldData.m_f0 * g0;
    grad += fieldData.m_f1 * g1;
    grad += fieldData.m_f2 * g2;

    for (pmp::Vertex vv : fieldData.m_g0)
    {
        J[vv] += g0 * fieldData.m_g0[vv];
    }
    for (pmp::Vertex vv : fieldData.m_g1)
    {
        J[vv] += g1 * fieldData.m_g1[vv];
    }
    for (pmp::Vertex vv : fieldData.m_g2)
    {
        J[vv] += g2 * fieldData.m_g2[vv];
    }

    const Eigen::Matrix3d cn = crossMatrix<double>(n);
    const Eigen::Matrix3d ce01 = crossMatrix<double>(e01);
    const Eigen::Matrix3d ce12 = crossMatrix<double>(e12);
    const Eigen::Matrix3d ce20 = crossMatrix<double>(e20);

    const Eigen::Matrix3d P = (Eigen::Matrix3d::Identity() - n * n.transpose()) / doubleArea;
    const Eigen::Matrix3d Jnv0 = P * ce12;
    const Eigen::Matrix3d Jnv1 = P * ce20;
    const Eigen::Matrix3d Jnv2 = P * ce01;

    J[v0] += fieldData.m_f0 * (-ce12 * Jnv0 - g0 * n.cross(e12).transpose()) / doubleArea;
    J[v1] += fieldData.m_f0 * (-cn - ce12 * Jnv1 - g0 * n.cross(e20).transpose()) / doubleArea;
    J[v2] += fieldData.m_f0 * (cn - ce12 * Jnv2 - g0 * n.cross(e01).transpose()) / doubleArea;

    J[v0] += fieldData.m_f1 * (cn - ce20 * Jnv0 - g1 * n.cross(e12).transpose()) / doubleArea;
    J[v1] += fieldData.m_f1 * (-ce20 * Jnv1 - g1 * n.cross(e20).transpose()) / doubleArea;
    J[v2] += fieldData.m_f1 * (-cn - ce20 * Jnv2 - g1 * n.cross(e01).transpose()) / doubleArea;

    J[v0] += fieldData.m_f2 * (-cn - ce01 * Jnv0 - g2 * n.cross(e12).transpose()) / doubleArea;
    J[v1] += fieldData.m_f2 * (cn - ce01 * Jnv1 - g2 * n.cross(e20).transpose()) / doubleArea;
    J[v2] += fieldData.m_f2 * (-ce01 * Jnv2 - g2 * n.cross(e01).transpose()) / doubleArea;

    return grad;
}
