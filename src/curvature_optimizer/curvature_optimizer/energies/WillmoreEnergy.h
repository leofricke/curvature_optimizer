#pragma once

#include "opt/GaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;

struct WillmoreEnergy : public GaussNewtonEnergyModule
{
    explicit WillmoreEnergy(CurvatureOptimizerWorkspace& ws);

    size_t numElements() const override;
    size_t numResiduals(size_t element) const override;
    void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const override;
    void initialize() override;

    CurvatureOptimizerWorkspace& m_ws;
};
