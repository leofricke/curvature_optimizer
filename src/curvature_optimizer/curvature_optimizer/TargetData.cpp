#include "TargetData.h"

#include "CurvatureOptimizerWorkspace.h"

TargetMetric createTargetMetric(const CurvatureOptimizerWorkspace& ws, pmp::Face f)
{
    TargetMetric metric;
    metric.m_face = f;
    metric.m_area = ws.m_currentFaceAreas[f];

    metric.m_subdivision.updateCenterWeights(ws, f);
    metric.m_subdivision.initVirtualTriangles(ws, f, metric.m_triangles);

    return metric;
}
