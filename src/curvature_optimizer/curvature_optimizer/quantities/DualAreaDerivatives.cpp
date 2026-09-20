#include "DualAreaDerivatives.h"

#include "DiffRules.h"
#include "../CurvatureOptimizerWorkspace.h"
#include "../VirtualSubdivision.h"

double computeDualArea(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, DualAreaMode mode)
{
    switch (mode)
    {
    case DualAreaMode::Barycentric:
        return computeBarycentricArea(ws, v);
    case DualAreaMode::Voronoi:
        return computeVoronoiArea(ws, v, false);
    case DualAreaMode::VoronoiMixed:
        return computeVoronoiArea(ws, v, true);
    };
    throw std::runtime_error("Invalid DualAreaMode");
}

double computeDualAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g, DualAreaMode mode)
{
    switch (mode)
    {
    case DualAreaMode::Barycentric:
        return computeBarycentricAreaGradient(ws, v, g);
    case DualAreaMode::Voronoi:
        return computeVoronoiAreaGradient(ws, v, g, false);
    case DualAreaMode::VoronoiMixed:
        return computeVoronoiAreaGradient(ws, v, g, true);
    };
    throw std::runtime_error("Invalid DualAreaMode");
}

double computeBarycentricArea(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v)
{
    double A = 0.0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(h))
        {
            const pmp::Face f = ws.m_mesh.face(h);
            const double val = static_cast<double>(ws.m_mesh.valence(f));

            A += ws.m_currentFaceAreas[f] / val;
        }
    }

    return A;
}

double computeBarycentricAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g)
{
    double A = 0.0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(h))
        {
            const pmp::Face f = ws.m_mesh.face(h);
            const double val = static_cast<double>(ws.m_mesh.valence(f));

            A += ws.m_currentFaceAreas[f] / val;

            const LocalGradient<Eigen::RowVector3d, 4>& ga = ws.m_currentFaceAreaGradients[f];

            g.multiplyAdd(1.0 / val, ga);
        }
    }

    return A;
}

double computeVoronoiArea(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, bool mixed)
{
    double A = 0.0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(h))
        {
            const pmp::Face f = ws.m_mesh.face(h);
            if (ws.m_mesh.valence(f) == 3)
            {
                const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
                const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
                const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
                const Eigen::Vector3d p0 = ws.position(v0);
                const Eigen::Vector3d p1 = ws.position(v1);
                const Eigen::Vector3d p2 = ws.position(v2);
                const Eigen::Vector3d e01 = p1 - p0;
                const Eigen::Vector3d e12 = p2 - p1;
                const Eigen::Vector3d e20 = p0 - p2;
                const double doubleArea = e01.cross(-e20).norm();
                const double d0 = e01.dot(-e20);
                const double d1 = e12.dot(-e01);
                const double d2 = e20.dot(-e12);
                const double le01 = e01.squaredNorm();
                const double le20 = e20.squaredNorm();

                if (mixed && d0 < 0.0)
                {
                    A += 0.25 * doubleArea;
                }
                else if (mixed && (d1 < 0.0 || d2 < 0.0))
                {
                    A += 0.125 * doubleArea;
                }
                else
                {
                    const double cot1 = d1 / doubleArea;
                    const double cot2 = d2 / doubleArea;
                    A += 0.125 * (le01 * cot2 + le20 * cot1);
                }
            }
            else
            {
                VirtualSubdivision subdiv;
                subdiv.updateCenterWeights(ws, f);
                const Eigen::Vector3d center = subdiv.getCenterPoint(ws, f);
                const double weight = subdiv.getVertexWeight(ws, f, v);

                for (pmp::Halfedge hh : ws.m_mesh.halfedges(f))
                {
                    const Eigen::Vector3d p0 = ws.position(ws.m_mesh.from_vertex(hh));
                    const Eigen::Vector3d p1 = ws.position(ws.m_mesh.to_vertex(hh));
                    const Eigen::Vector3d e01 = p1 - p0;
                    const Eigen::Vector3d e12 = center - p1;
                    const Eigen::Vector3d e20 = p0 - center;
                    const double doubleArea = e01.cross(-e20).norm();
                    const double d0 = e01.dot(-e20);
                    const double d1 = e12.dot(-e01);
                    const double d2 = e20.dot(-e12);
                    const double cot0 = d0 / doubleArea;
                    const double cot1 = d1 / doubleArea;
                    const double cot2 = d2 / doubleArea;
                    const double le01 = e01.squaredNorm();
                    const double le12 = e12.squaredNorm();
                    const double le20 = e20.squaredNorm();

                    // outer voronoi cell
                    if (ws.m_mesh.from_vertex(hh) == v)
                    {
                        if (mixed && d0 < 0.0)
                        {
                            A += 0.25 * doubleArea;
                        }
                        else if (mixed && (d1 < 0.0 || d2 < 0.0))
                        {
                            A += 0.125 * doubleArea;
                        }
                        else
                        {
                            A += 0.125 * (le01 * cot2 + le20 * cot1);
                        }
                    }
                    else if (ws.m_mesh.to_vertex(hh) == v)
                    {
                        if (mixed && d1 < 0.0)
                        {
                            A += 0.25 * doubleArea;
                        }
                        else if (mixed && (d0 < 0.0 || d2 < 0.0))
                        {
                            A += 0.125 * doubleArea;
                        }
                        else
                        {
                            A += 0.125 * (le12 * cot0 + le01 * cot2);
                        }
                    }

                    // inner voronoi cell
                    if (mixed && d2 < 0.0)
                    {
                        A += weight * 0.25 * doubleArea;
                    }
                    else if (mixed && (d0 < 0.0 || d1 < 0.0))
                    {
                        A += weight * 0.125 * doubleArea;
                    }
                    else
                    {
                        A += weight * 0.125 * (le20 * cot1 + le12 * cot0);
                    }
                }
            }
        }
    }

    return A;
}

double computeVoronoiAreaGradient(const CurvatureOptimizerWorkspace& ws, pmp::Vertex v, LocalGradient<Eigen::RowVector3d, 8>& g, bool mixed)
{
    double A = 0.0;
    for (pmp::Halfedge h : ws.m_mesh.halfedges(v))
    {
        if (!ws.m_mesh.is_boundary(h))
        {
            const pmp::Face f = ws.m_mesh.face(h);
            if (ws.m_mesh.valence(f) == 3)
            {
                const pmp::Vertex v0 = ws.m_mesh.from_vertex(h);
                const pmp::Vertex v1 = ws.m_mesh.to_vertex(h);
                const pmp::Vertex v2 = ws.m_mesh.to_vertex(ws.m_mesh.next_halfedge(h));
                const Eigen::Vector3d p0 = ws.position(v0);
                const Eigen::Vector3d p1 = ws.position(v1);
                const Eigen::Vector3d p2 = ws.position(v2);

                const Eigen::Vector3d e01 = p1 - p0;
                const Eigen::Vector3d e12 = p2 - p1;
                const Eigen::Vector3d e20 = p0 - p2;
                const double doubleArea = e01.cross(-e20).norm();
                const double d0 = e01.dot(-e20);
                const double d1 = e12.dot(-e01);
                const double d2 = e20.dot(-e12);
                const LocalGradient<Eigen::RowVector3d, 4>& ga = ws.m_currentFaceAreaGradients[ws.m_mesh.face(h)];

                if (mixed && d0 < 0.0)
                {
                    for (pmp::Vertex vv : ga)
                    {
                        g[vv] += 0.5 * ga[vv];
                    }
                    A += 0.25 * doubleArea;
                }
                else if (mixed && (d1 < 0.0 || d2 < 0.0))
                {
                    for (pmp::Vertex vv : ga)
                    {
                        g[vv] += 0.25 * ga[vv];
                    }
                    A += 0.125 * doubleArea;
                }
                else
                {
                    const double cot1 = d1 / doubleArea;
                    const double cot2 = d2 / doubleArea;
                    const double le01 = e01.squaredNorm();
                    const double le20 = e20.squaredNorm();

                    const Eigen::RowVector3d dcot1v0 = (e12.transpose() + d1 * -2.0 * ga[v0] / doubleArea) / doubleArea;
                    const Eigen::RowVector3d dcot1v2 = (-e01.transpose() + d1 * -2.0 * ga[v2] / doubleArea) / doubleArea;
                    const Eigen::RowVector3d dcot1v1 = -dcot1v0 - dcot1v2;
                    const Eigen::RowVector3d dcot2v0 = (-e12.transpose() + d2 * -2.0 * ga[v0] / doubleArea) / doubleArea;
                    const Eigen::RowVector3d dcot2v1 = (e20.transpose() + d2 * -2.0 * ga[v1] / doubleArea) / doubleArea;
                    const Eigen::RowVector3d dcot2v2 = -dcot2v0 - dcot2v1;

                    g[v0] += 0.125 * (-2.0 * e01.transpose() * cot2 + 2.0 * e20.transpose() * cot1 + le20 * dcot1v0 + le01 * dcot2v0);
                    g[v1] += 0.125 * (2.0 * e01.transpose() * cot2 + le20 * dcot1v1 + le01 * dcot2v1);
                    g[v2] += 0.125 * (-2.0 * e20.transpose() * cot1 + le20 * dcot1v2 + le01 * dcot2v2);

                    A += 0.125 * (le01 * cot2 + le20 * cot1);
                }
            }
            else
            {
                VirtualSubdivision subdiv;
                subdiv.updateCenterWeights(ws, f);
                const Eigen::Vector3d center = subdiv.getCenterPoint(ws, f);
                const double weight = subdiv.getVertexWeight(ws, f, v);

                for (pmp::Halfedge hh : ws.m_mesh.halfedges(f))
                {
                    const pmp::Vertex v0 = ws.m_mesh.from_vertex(hh);
                    const pmp::Vertex v1 = ws.m_mesh.to_vertex(hh);
                    const Eigen::Vector3d p0 = ws.position(v0);
                    const Eigen::Vector3d p1 = ws.position(v1);
                    const Eigen::Vector3d e01 = p1 - p0;
                    const Eigen::Vector3d e12 = center - p1;
                    const Eigen::Vector3d e20 = p0 - center;
                    Eigen::Vector3d n = e01.cross(-e20);
                    const double doubleArea = n.norm();
                    n.normalize();
                    const double d0 = e01.dot(-e20);
                    const double d1 = e12.dot(-e01);
                    const double d2 = e20.dot(-e12);
                    const double cot0 = d0 / doubleArea;
                    const double cot1 = d1 / doubleArea;
                    const double cot2 = d2 / doubleArea;
                    const double le01 = e01.squaredNorm();
                    const double le12 = e12.squaredNorm();
                    const double le20 = e20.squaredNorm();

                    LocalGradient<Eigen::RowVector3d, 8> ga;
                    ga[v0] += n.cross(e12);
                    ga[v1] += n.cross(e20);
                    {
                        size_t hidx = 0;
                        for (pmp::Halfedge hhh : ws.m_mesh.halfedges(f))
                        {
                            ga[ws.m_mesh.from_vertex(hhh)] += subdiv.m_centerWeights[hidx++] * n.cross(e01);
                        }
                    }

                    // cot derivs
                    LocalGradient<Eigen::RowVector3d, 8> gcot0;
                    LocalGradient<Eigen::RowVector3d, 8> gcot1;
                    LocalGradient<Eigen::RowVector3d, 8> gcot2;
                    gcot0[v0] += e20 - e01;
                    gcot0[v1] += -e20;
                    gcot1[v0] += e12;
                    gcot1[v1] += -e12 + e01;
                    gcot2[v0] += -e12;
                    gcot2[v1] += e20;
                    {
                        size_t hidx = 0;
                        for (pmp::Halfedge hhh : ws.m_mesh.halfedges(f))
                        {
                            gcot0[ws.m_mesh.from_vertex(hhh)] += subdiv.m_centerWeights[hidx] * e01;
                            gcot1[ws.m_mesh.from_vertex(hhh)] += subdiv.m_centerWeights[hidx] * -e01;
                            gcot2[ws.m_mesh.from_vertex(hhh)] += subdiv.m_centerWeights[hidx] * (e12 - e20);
                            hidx++;
                        }
                    }
                    applyQuotientRule1Inplace(d0, doubleArea, ga, gcot0);
                    applyQuotientRule1Inplace(d1, doubleArea, ga, gcot1);
                    applyQuotientRule1Inplace(d2, doubleArea, ga, gcot2);

                    // lengths derivs
                    LocalGradient<Eigen::RowVector3d, 8> gle01;
                    LocalGradient<Eigen::RowVector3d, 8> gle12;
                    LocalGradient<Eigen::RowVector3d, 8> gle20;
                    gle01[v0] += -2.0 * e01;
                    gle01[v1] += 2.0 * e01;
                    gle12[v1] += -2.0 * e12;
                    gle20[v0] += 2.0 * e20;
                    {
                        size_t hidx = 0;
                        for (pmp::Halfedge hhh : ws.m_mesh.halfedges(f))
                        {
                            gle12[ws.m_mesh.from_vertex(hhh)] += subdiv.m_centerWeights[hidx] * 2.0 * e12;
                            gle20[ws.m_mesh.from_vertex(hhh)] += subdiv.m_centerWeights[hidx] * -2.0 * e20;
                            hidx++;
                        }
                    }

                    // outer voronoi cell
                    if (ws.m_mesh.from_vertex(hh) == v)
                    {
                        if (mixed && d0 < 0.0)
                        {
                            for (pmp::Vertex vv : ga)
                            {
                                g[vv] += 0.25 * ga[vv];
                            }
                            A += 0.25 * doubleArea;
                        }
                        else if (mixed && (d1 < 0.0 || d2 < 0.0))
                        {
                            for (pmp::Vertex vv : ga)
                            {
                                g[vv] += 0.125 * ga[vv];
                            }
                            A += 0.125 * doubleArea;
                        }
                        else
                        {
                            for (pmp::Vertex vv : gcot2)
                            {
                                g[vv] += 0.125 * le01 * gcot2[vv];
                            }
                            for (pmp::Vertex vv : gcot1)
                            {
                                g[vv] += 0.125 * le20 * gcot1[vv];
                            }
                            for (pmp::Vertex vv : gle01)
                            {
                                g[vv] += 0.125 * gle01[vv] * cot2;
                            }
                            for (pmp::Vertex vv : gle20)
                            {
                                g[vv] += 0.125 * gle20[vv] * cot1;
                            }
                            A += 0.125 * (le01 * cot2 + le20 * cot1);
                        }
                    }
                    else if (ws.m_mesh.to_vertex(hh) == v)
                    {
                        if (mixed && d1 < 0.0)
                        {
                            for (pmp::Vertex vv : ga)
                            {
                                g[vv] += 0.25 * ga[vv];
                            }
                            A += 0.25 * doubleArea;
                        }
                        else if (mixed && (d0 < 0.0 || d2 < 0.0))
                        {
                            for (pmp::Vertex vv : ga)
                            {
                                g[vv] += 0.125 * ga[vv];
                            }
                            A += 0.125 * doubleArea;
                        }
                        else
                        {
                            for (pmp::Vertex vv : gcot0)
                            {
                                g[vv] += 0.125 * le12 * gcot0[vv];
                            }
                            for (pmp::Vertex vv : gcot2)
                            {
                                g[vv] += 0.125 * le01 * gcot2[vv];
                            }
                            for (pmp::Vertex vv : gle12)
                            {
                                g[vv] += 0.125 * gle12[vv] * cot0;
                            }
                            for (pmp::Vertex vv : gle01)
                            {
                                g[vv] += 0.125 * gle01[vv] * cot2;
                            }
                            A += 0.125 * (le12 * cot0 + le01 * cot2);
                        }
                    }

                    // inner voronoi cell
                    if (mixed && d2 < 0.0)
                    {
                        for (pmp::Vertex vv : ga)
                        {
                            g[vv] += weight * 0.25 * ga[vv];
                        }
                        A += weight * 0.25 * doubleArea;
                    }
                    else if (mixed && (d0 < 0.0 || d1 < 0.0))
                    {
                        for (pmp::Vertex vv : ga)
                        {
                            g[vv] += weight * 0.125 * ga[vv];
                        }
                        A += weight * 0.125 * doubleArea;
                    }
                    else
                    {
                        for (pmp::Vertex vv : gcot1)
                        {
                            g[vv] += weight * 0.125 * le20 * gcot1[vv];
                        }
                        for (pmp::Vertex vv : gcot0)
                        {
                            g[vv] += weight * 0.125 * le12 * gcot0[vv];
                        }
                        for (pmp::Vertex vv : gle20)
                        {
                            g[vv] += weight * 0.125 * gle20[vv] * cot1;
                        }
                        for (pmp::Vertex vv : gle12)
                        {
                            g[vv] += weight * 0.125 * gle12[vv] * cot0;
                        }
                        A += weight * 0.125 * (le20 * cot1 + le12 * cot0);
                    }
                }
            }
        }
    }

    return A;
}
