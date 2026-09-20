#include "CurvatureTensorVertexCurvature.h"

#include "CurvatureOptimizerWorkspace.h"
#include "quantities/CommonDerivatives.h"
#include "utils/KoppEigendecomposition.h"
#include "Tolerances.h"

CurvatureTensorVertexCurvature::CurvatureTensorVertexCurvature()
{
}

Eigen::Matrix3d CurvatureTensorVertexCurvature::computeTensor(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v) const
{
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();

    const double A = ws.m_currentVertexDualAreas[v];

    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        const pmp::Edge e = ws.m_mesh.edge(h);
        if (!ws.m_mesh.is_boundary(e))
        {
            const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
            const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
            const Eigen::Vector3d p0 = ws.position(v0);
            const Eigen::Vector3d p1 = ws.position(v1);

            const Eigen::Vector3d ev = p1 - p0;
            const double l = ev.norm();
            const Eigen::Vector3d ev0 = ev.normalized();

            const double beta = ws.m_currentDihedralAngles[e];
            const Eigen::Matrix3d ee = ev0 * ev0.transpose();

            tensor += beta * 0.5 * l * ee;
        }
    }
    tensor /= A;

    return tensor;
}

Eigen::Matrix3d CurvatureTensorVertexCurvature::computeTensorDeriv(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::Matrix<double, 3, 9>>& dT) const
{
    Eigen::Matrix3d tensor = Eigen::Matrix3d::Zero();

    const double A = ws.m_currentVertexDualAreas[v];
    const LocalGradient<Eigen::RowVector3d, 8>& ga = ws.m_currentVertexDualAreaGradients[v];

    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        const pmp::Edge e = ws.m_mesh.edge(h);
        if (!ws.m_mesh.is_boundary(e))
        {
            const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
            const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
            const Eigen::Vector3d p0 = ws.position(v0);
            const Eigen::Vector3d p1 = ws.position(v1);

            const Eigen::Vector3d ev = p1 - p0;
            const double l = ev.norm();
            const Eigen::Vector3d ev0 = ev.normalized();

            const double beta = ws.m_currentDihedralAngles[e];
            const Eigen::Matrix3d ee = ev0 * ev0.transpose();

            tensor += beta * 0.5 * l * ee;

            // dihedral angle term
            const LocalGradient<Eigen::RowVector3d>& gb = ws.m_currentDihedralAngleGradients[e];
            for (pmp::Vertex vv : gb)
            {
                for (size_t i = 0; i < 3; ++i)
                {
                    dT[vv].block<3, 3>(0, 3 * i) += gb[vv][i] * 0.5 * l * ee;
                }
            }

            // length term
            for (size_t i = 0; i < 3; ++i)
            {
                dT[v0].block<3, 3>(0, 3 * i) += beta * 0.5 * -ev0[i] * ee;
                dT[v1].block<3, 3>(0, 3 * i) += beta * 0.5 * ev0[i] * ee;
            }

            // normalized outer product term
            std::array<Eigen::Matrix3d, 3> gee;
            computeNormalizedOuterProductJacobian(ev, gee);

            for (size_t i = 0; i < 3; ++i)
            {
                dT[v0].block<3, 3>(0, 3 * i) += beta * 0.5 * l * -gee[i];
                dT[v1].block<3, 3>(0, 3 * i) += beta * 0.5 * l * gee[i];
            }
        }
    }

    tensor /= A;
    dT /= A;
    for (pmp::Vertex vv : ga)
    {
        for (size_t i = 0; i < 3; ++i)
        {
            dT[vv].block<3, 3>(0, 3 * i) -= tensor * ga[vv][i] / A;
        }
    }

    return tensor;
}

Eigen::Vector2d CurvatureTensorVertexCurvature::computePrincipalCurvatures(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    const Eigen::Matrix3d tensor = computeTensor(ws, v);
    updateCurvaturesAndDirections(tensor);

    return { m_curv.m_kmax, m_curv.m_kmin };
}

Eigen::Vector2d CurvatureTensorVertexCurvature::computePrincipalCurvaturesGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& gkmax,
                                                                                   LocalGradient<Eigen::RowVector3d, 8>& gkmin)
{
    LocalGradient<Eigen::Matrix<double, 3, 9>> dT;
    const Eigen::Matrix3d tensor = computeTensorDeriv(ws, v, dT);
    updateCurvaturesAndDirections(tensor);

    for (pmp::Vertex vv : dT)
    {
        for (size_t i = 0; i < 3; ++i)
        {
            const Eigen::Matrix3d& dTvv = dT[vv].block<3, 3>(0, 3 * i);
            const double trace = dTvv.trace();

            Eigen::Vector2d dc;
            if (std::abs(m_curv.m_kmax) < EpsCurv && std::abs(m_curv.m_kmin) >= EpsCurv) // cylindric point
            {
                dc[1] = m_curv.m_minDir.dot(dTvv * m_curv.m_minDir);
                dc[0] = trace - dc[1];
            }
            else if (std::abs(m_curv.m_kmin) < EpsCurv && std::abs(m_curv.m_kmax) >= EpsCurv) // cylindric point
            {
                dc[0] = m_curv.m_maxDir.dot(dTvv * m_curv.m_maxDir);
                dc[1] = trace - dc[0];
            }
            else if (std::abs(m_curv.m_kmax - m_curv.m_kmin) < EpsCurv) // umbilic point
            {
                dc[0] = dc[1] = 0.5 * trace;
            }
            else // distinct curvatures
            {
                dc[0] = m_curv.m_maxDir.dot(dTvv * m_curv.m_maxDir);
                dc[1] = m_curv.m_minDir.dot(dTvv * m_curv.m_minDir);
            }

            gkmax[vv][i] += dc[0];
            gkmin[vv][i] += dc[1];
        }
    }

    return { m_curv.m_kmax, m_curv.m_kmin };
}

void CurvatureTensorVertexCurvature::updateCurvaturesAndDirections(const Eigen::Matrix3d& tensor)
{
    Eigen::Vector3d w;
    Eigen::Matrix3d Q;
    eigendecomposition3(tensor, w, Q);

    const double eval0 = w[2];
    const double eval1 = w[1];
    const double eval2 = w[0];
    Eigen::Ref<const Eigen::Vector3d> evec0 = Q.col(2);
    Eigen::Ref<const Eigen::Vector3d> evec1 = Q.col(1);
    Eigen::Ref<const Eigen::Vector3d> evec2 = Q.col(0);

    const double a0 = std::abs(eval0);
    const double a1 = std::abs(eval1);
    const double a2 = std::abs(eval2);
    if (a0 < a1)
    {
        if (a0 < a2)
        {
            m_curv.m_kmax = eval1;
            m_curv.m_kmin = eval2;
            m_curv.m_maxDir = evec1;
            m_curv.m_minDir = evec2;
        }
        else
        {
            m_curv.m_kmax = eval0;
            m_curv.m_kmin = eval1;
            m_curv.m_maxDir = evec0;
            m_curv.m_minDir = evec1;
        }
    }
    else
    {
        if (a1 < a2)
        {
            m_curv.m_kmax = eval0;
            m_curv.m_kmin = eval2;
            m_curv.m_maxDir = evec0;
            m_curv.m_minDir = evec2;
        }
        else
        {
            m_curv.m_kmax = eval0;
            m_curv.m_kmin = eval1;
            m_curv.m_maxDir = evec0;
            m_curv.m_minDir = evec1;
        }
    }
}
