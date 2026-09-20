#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;
struct TargetMetric;

struct ConformalEnergy : public ModularGaussNewtonEnergyBlock
{
    ConformalEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetMetric>& targetMetrics);

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    const CurvatureOptimizerWorkspace& m_ws;
    const std::vector<TargetMetric>& m_targetMetrics;
};
