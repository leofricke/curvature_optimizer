#pragma once

#include "../CurvatureOptimizerWorkspace.h"
#include "opt/GaussNewtonEnergyFunction.h"

/**
* Inserts a per-vertex gradient (e.g. a LocalGradient) into a Jacobian row, skipping locked vertices.
*/
template <typename Gradient>
void insertVertexGradient(JacobianRow& row, const CurvatureOptimizerWorkspace& ws, const Gradient& g, double scale = 1.0)
{
    for (pmp::Vertex v : g)
    {
        const int dofIdx = ws.m_dofs.dofIdx(v);
        if (dofIdx >= 0)
        {
            row.insertSegment(dofIdx * 3, scale * g[v]);
        }
    }
}
