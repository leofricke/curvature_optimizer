#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;
struct TargetPosition;

struct PositionEnergy : public ModularGaussNewtonEnergyBlock
{
    PositionEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetPosition>& targetPositions);

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    const CurvatureOptimizerWorkspace& m_ws;
    const std::vector<TargetPosition>& m_targetPositions;
};
