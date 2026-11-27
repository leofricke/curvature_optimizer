#pragma once

#include "opt/GaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;
struct TargetVertexCurvature;

struct CurvatureEnergy : public GaussNewtonEnergyModule
{
    CurvatureEnergy(CurvatureOptimizerWorkspace& ws, const std::vector<TargetVertexCurvature>& targetVertexCurvatures);

    size_t numElements() const override;
    size_t numResiduals(size_t element) const override;
    void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const override;
    void initialize() override;

    CurvatureOptimizerWorkspace& m_ws;
    const std::vector<TargetVertexCurvature>& m_targetVertexCurvatures;
};
