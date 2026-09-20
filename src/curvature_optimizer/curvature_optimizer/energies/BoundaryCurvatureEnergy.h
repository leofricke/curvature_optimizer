#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;
struct TargetBoundaryCurvature;

struct BoundaryCurvatureEnergy : public ModularGaussNewtonEnergyBlock
{
    BoundaryCurvatureEnergy(const CurvatureOptimizerWorkspace& ws, const std::vector<TargetBoundaryCurvature>& targetBoundaryCurvatures);

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    const CurvatureOptimizerWorkspace& m_ws;
    const std::vector<TargetBoundaryCurvature>& m_targetBoundaryCurvatures;
};
