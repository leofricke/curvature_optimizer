#include "TargetData.h"

TargetMetric createTargetMetric(const pmp::SurfaceMesh& mesh, pmp::Face f)
{
    TargetMetric metric;
    metric.m_face = f;

    Eigen::Vector3d areaVector = Eigen::Vector3d::Zero();
    for (pmp::Halfedge h : mesh.halfedges(f))
    {
        const Eigen::Vector3d p0 = static_cast<Eigen::Vector3d>(mesh.position(mesh.from_vertex(h)));
        const Eigen::Vector3d p1 = static_cast<Eigen::Vector3d>(mesh.position(mesh.to_vertex(h)));
        areaVector += p0.cross(p1);
    }
    metric.m_area = 0.5 * areaVector.norm();

    metric.m_subdivision.updateCenterWeights(mesh, f);
    metric.m_subdivision.initVirtualTriangles(mesh, f, metric.m_triangles);

    return metric;
}
