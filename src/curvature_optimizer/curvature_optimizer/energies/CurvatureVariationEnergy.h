#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;

struct CurvatureVariationEnergy : public ModularGaussNewtonEnergyBlock
{
    CurvatureVariationEnergy(CurvatureOptimizerWorkspace& ws, double pNorm = 2.0);

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    CurvatureOptimizerWorkspace& m_ws;
    double m_pNorm;
};
