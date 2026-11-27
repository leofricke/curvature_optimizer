#pragma once

#include "opt/GaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;
struct TargetPosition;

struct PositionEnergy : public GaussNewtonEnergyModule
{
    PositionEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetPosition>& targetPositions);

    size_t numElements() const override;
    size_t numResiduals(size_t element) const override;
    void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const override;
    void initialize() override;

    const CurvatureOptimizerWorkspace& m_ws;
    const std::vector<TargetPosition>& m_targetPositions;
};
