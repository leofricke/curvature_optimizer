#pragma once

#include "opt/GaussNewtonEnergyFunction.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

struct PolygonPlanarityEnergy : public GaussNewtonEnergyModule
{
    PolygonPlanarityEnergy(const CurvatureOptimizerWorkspace& ws);

    size_t numElements() const override;
    size_t numResiduals(size_t element) const override;
    void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const override;
    void initialize() override;

    const CurvatureOptimizerWorkspace& m_ws;
    std::vector<pmp::Face> m_polygons;
};
