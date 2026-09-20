#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;

struct PolygonPlanarityEnergy : public ModularGaussNewtonEnergyBlock
{
    PolygonPlanarityEnergy(const CurvatureOptimizerWorkspace& ws);

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    const CurvatureOptimizerWorkspace& m_ws;
};
