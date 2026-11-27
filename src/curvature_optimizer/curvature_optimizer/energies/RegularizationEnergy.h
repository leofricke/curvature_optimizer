#pragma once

#include "opt/GaussNewtonEnergyFunction.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

struct RegularizationEnergy : public GaussNewtonEnergyModule
{
    struct InitialPosition
    {
        pmp::Vertex m_vertex;
        Eigen::Vector3d m_position;
    };

    RegularizationEnergy(const CurvatureOptimizerWorkspace& ws);

    void resetPositions();

    size_t numElements() const override;
    size_t numResiduals(size_t element) const override;
    void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const override;
    void initialize() override;

    const CurvatureOptimizerWorkspace& m_ws;
    std::vector<InitialPosition> m_initialPositions;
};
