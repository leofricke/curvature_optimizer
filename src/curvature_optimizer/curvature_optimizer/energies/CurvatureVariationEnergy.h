#pragma once

#include "opt/GaussNewtonEnergyFunction.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

struct CurvatureVariationEnergy : public GaussNewtonEnergyModule
{
    CurvatureVariationEnergy(CurvatureOptimizerWorkspace& ws, double pNorm = 2.0);

    size_t numElements() const override;
    size_t numResiduals(size_t element) const override;
    void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const override;

    CurvatureOptimizerWorkspace& m_ws;
    double m_pNorm;
    std::vector<pmp::Face> m_faces; // faces without boundary vertices
};
