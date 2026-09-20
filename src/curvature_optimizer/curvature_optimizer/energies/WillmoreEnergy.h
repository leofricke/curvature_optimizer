#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;

struct WillmoreEnergy : public ModularGaussNewtonEnergyBlock
{
    explicit WillmoreEnergy(CurvatureOptimizerWorkspace& ws);

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    CurvatureOptimizerWorkspace& m_ws;
};
