#include "DihedralAngleDerivatives.h"

#include "DiffRules.h"
#include "FaceAreaDerivatives.h"
#include "../CurvatureOptimizerWorkspace.h"
#include "utils/EigenMathUtil.h"

// notation from: Rasmus Tamstorf, Eitan Grinspun "Discrete bending forces and their Jacobians"

double computeDihedralAngle(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, bool robustAngle)
{
    if (ws.m_mesh.is_boundary(e))
    {
        return 0.0;
    }
    double beta;
    if (ws.m_mesh.valence(ws.m_mesh.face(e, 0)) == 3 && ws.m_mesh.valence(ws.m_mesh.face(e, 1)) == 3)
    {
        beta = computeTriangleDihedralAngle(ws, e);
    }
    else
    {
        beta = computePolygonDihedralAngle(ws, e);
    }

    if (robustAngle)
    {
        beta = 2.0 * std::tan(0.5 * beta);
    }
    return beta;
}

double computeDihedralAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, LocalGradient<Eigen::RowVector3d, 4>& g, bool robustAngle)
{
    if (ws.m_mesh.is_boundary(e))
    {
        return 0.0;
    }
    double beta;
    if (ws.m_mesh.valence(ws.m_mesh.face(e, 0)) == 3 && ws.m_mesh.valence(ws.m_mesh.face(e, 1)) == 3)
    {
        beta = computeTriangleDihedralAngleGradient(ws, e, g);
    }
    else
    {
        beta = computePolygonDihedralAngleGradient(ws, e, g);
    }

    if (robustAngle)
    {
        const auto d1 = [](double x)
        {
            return 1.0 / std::pow(std::cos(0.5 * x), 2);
        };
        applyChainRule1Inplace(beta, d1, g);
        beta = 2.0 * std::tan(0.5 * beta);
    }
    return beta;
}

double computeTriangleDihedralAngle(const CurvatureOptimizerWorkspace& ws, pmp::Edge e)
{
    if (ws.m_mesh.is_boundary(e))
    {
        return 0.0;
    }
    const pmp::Halfedge h = ws.m_mesh.halfedge(e, 0);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d e01 = p1 - p0;
    const double le01 = e01.norm();

    const Eigen::Vector3d& n0 = ws.m_currentFaceNormals[ws.m_mesh.face(h)];
    const Eigen::Vector3d& n1 = ws.m_currentFaceNormals[ws.m_mesh.face(ws.m_mesh.opposite_halfedge(h))];

    const double sin = std::clamp(e01.dot(n0.cross(n1)) / le01, -1.0, 1.0);
    const double cos = std::clamp(n0.dot(n1), -1.0, 1.0);

    return std::atan2(sin, cos);
}

double computeTriangleDihedralAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, LocalGradient<Eigen::RowVector3d, 4>& g)
{
    if (ws.m_mesh.is_boundary(e))
    {
        return 0.0;
    }
    const pmp::Halfedge h = ws.m_mesh.halfedge(e, 0);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
    const pmp::Vertex v3 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(ws.m_mesh.opposite_halfedge(h)));
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d p2 = ws.position(v2);
    const Eigen::Vector3d p3 = ws.position(v3);
    const Eigen::Vector3d e01 = p1 - p0;
    const Eigen::Vector3d e02 = p2 - p0;
    const Eigen::Vector3d e03 = p3 - p0;
    const Eigen::Vector3d e12 = p2 - p1;
    const Eigen::Vector3d e13 = p3 - p1;
    const double le01 = e01.norm();
    const double le02 = e02.norm();
    const double le03 = e03.norm();
    const double le12 = e12.norm();
    const double le13 = e13.norm();

    const Eigen::Vector3d& n0 = ws.m_currentFaceNormals[ws.m_mesh.face(h)];
    const Eigen::Vector3d& n1 = ws.m_currentFaceNormals[ws.m_mesh.face(ws.m_mesh.opposite_halfedge(h))];
    const double doubleArea0 = 2.0 * ws.m_currentFaceAreas[ws.m_mesh.face(h)];
    const double doubleArea1 = 2.0 * ws.m_currentFaceAreas[ws.m_mesh.face(ws.m_mesh.opposite_halfedge(h))];

    // inverse height
    std::array<double, 3> h0, h1;
    h0[0] = le12 / doubleArea0;
    h0[1] = le02 / doubleArea0;
    h0[2] = le01 / doubleArea0;
    h1[0] = le13 / doubleArea1;
    h1[1] = le03 / doubleArea1;
    h1[2] = le01 / doubleArea1;

    // cosines
    std::array<double, 2> c0, c1;
    c0[1] = e12.dot(-e01) / (le12 * le01);
    c0[0] = e02.dot(e01) / (le02 * le01);
    c1[1] = e13.dot(-e01) / (le13 * le01);
    c1[0] = e03.dot(e01) / (le03 * le01);

    g[v2] = -n0 * h0[2];
    g[v3] = -n1 * h1[2];
    g[v0] = n0 * h0[0] * c0[1] + n1 * h1[0] * c1[1];
    g[v1] = n0 * h0[1] * c0[0] + n1 * h1[1] * c1[0];

    const double sin = std::clamp(e01.dot(n0.cross(n1)) / le01, -1.0, 1.0);
    const double cos = std::clamp(n0.dot(n1), -1.0, 1.0);

    return std::atan2(sin, cos);
}

double computePolygonDihedralAngle(const CurvatureOptimizerWorkspace& ws, pmp::Edge e)
{
    if (ws.m_mesh.is_boundary(e))
    {
        return 0.0;
    }
    const pmp::Halfedge h = ws.m_mesh.halfedge(e, 0);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d e01 = p1 - p0;
    const double le01 = e01.norm();

    Eigen::Vector3d n0 = ws.m_currentFaceNormals[ws.m_mesh.face(h)];
    Eigen::Vector3d n1 = ws.m_currentFaceNormals[ws.m_mesh.face(ws.m_mesh.opposite_halfedge(h))];

    const Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - e01 * e01.transpose() / (le01 * le01);
    n0 = (P * n0).normalized();
    n1 = (P * n1).normalized();

    const double sin = std::clamp(e01.dot(n0.cross(n1)) / le01, -1.0, 1.0);
    const double cos = std::clamp(n0.dot(n1), -1.0, 1.0);

    return std::atan2(sin, cos);
}

double computePolygonDihedralAngleGradient(const CurvatureOptimizerWorkspace& ws, pmp::Edge e, LocalGradient<Eigen::RowVector3d, 4>& g)
{
    if (ws.m_mesh.is_boundary(e))
    {
        return 0.0;
    }
    const pmp::Halfedge h = ws.m_mesh.halfedge(e, 0);
    const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
    const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
    const Eigen::Vector3d p0 = ws.position(v0);
    const Eigen::Vector3d p1 = ws.position(v1);
    const Eigen::Vector3d e01 = p1 - p0;
    const double le01 = e01.norm();

    const Eigen::Matrix3d P = Eigen::Matrix3d::Identity() - e01 * e01.transpose() / (le01 * le01);
    LocalGradient<Eigen::Matrix3d> Jn0;
    LocalGradient<Eigen::Matrix3d> Jn1;
    Eigen::Vector3d n0 = computePolygonNormalJacobian(ws, ws.m_mesh.face(h), Jn0);
    Eigen::Vector3d n1 = computePolygonNormalJacobian(ws, ws.m_mesh.face(ws.m_mesh.opposite_halfedge(h)), Jn1);
    for (pmp::Vertex vv : Jn0)
    {
        Jn0[vv] = P * Jn0[vv];
    }
    const Eigen::Matrix3d de01 = (Eigen::Matrix3d::Identity() - e01 * e01.transpose() / (le01 * le01)) / le01;
    const Eigen::Matrix3d dPn0 = -(e01.dot(n0) * Eigen::Matrix3d::Identity() + e01 * n0.transpose()) / le01;
    Jn0[v0] -= dPn0 * de01;
    Jn0[v1] += dPn0 * de01;
    for (pmp::Vertex vv : Jn1)
    {
        Jn1[vv] = P * Jn1[vv];
    }
    const Eigen::Matrix3d dPn1 = -(e01.dot(n1) * Eigen::Matrix3d::Identity() + e01 * n1.transpose()) / le01;
    Jn1[v0] -= dPn1 * de01;
    Jn1[v1] += dPn1 * de01;

    n0 = P * n0;
    n1 = P * n1;

    const double ln0 = n0.norm();
    const double ln1 = n1.norm();
    const Eigen::Matrix3d Pn0 = (Eigen::Matrix3d::Identity() - n0 * n0.transpose() / (ln0 * ln0)) / ln0;
    const Eigen::Matrix3d Pn1 = (Eigen::Matrix3d::Identity() - n1 * n1.transpose() / (ln1 * ln1)) / ln1;
    for (pmp::Vertex vv : Jn0)
    {
        Jn0[vv] = Pn0 * Jn0[vv];
    }
    for (pmp::Vertex vv : Jn1)
    {
        Jn1[vv] = Pn1 * Jn1[vv];
    }
    n0.normalize();
    n1.normalize();

    const double sin = std::clamp(e01.dot(n0.cross(n1)) / le01, -1.0, 1.0);
    const double cos = std::clamp(n0.dot(n1), -1.0, 1.0);

    const Eigen::Matrix3d cn0 = crossMatrix<double>(n0);
    const Eigen::Matrix3d cn1 = crossMatrix<double>(n1);

    LocalGradient<Eigen::RowVector3d> gsin;
    LocalGradient<Eigen::RowVector3d> gcos;
    for (pmp::Vertex v : Jn0)
    {
        gsin[v] += (-cn1 * Jn0[v]).transpose() * e01 / le01;
        gcos[v] += Jn0[v].transpose() * n1;
    }
    for (pmp::Vertex v : Jn1)
    {
        gsin[v] += (cn0 * Jn1[v]).transpose() * e01 / le01;
        gcos[v] += Jn1[v].transpose() * n0;
    }

    const Eigen::RowVector3d dn = ((Eigen::Matrix3d::Identity() - e01 * e01.transpose() / (le01 * le01)) / le01) * n0.cross(n1);
    gsin[v1] += dn;
    gsin[v0] -= dn;

    for (pmp::Vertex v : gsin)
    {
        g[v] = gsin[v] * cos - sin * gcos[v];
    }

    return std::atan2(sin, cos);
}
