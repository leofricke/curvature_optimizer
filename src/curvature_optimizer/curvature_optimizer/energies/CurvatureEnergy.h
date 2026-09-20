#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"

struct CurvatureOptimizerWorkspace;
struct TargetVertexCurvature;

struct CurvatureEnergy : public ModularGaussNewtonEnergyBlock
{
    CurvatureEnergy(CurvatureOptimizerWorkspace& ws, const std::vector<TargetVertexCurvature>& targetVertexCurvatures);

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    CurvatureOptimizerWorkspace& m_ws;
    const std::vector<TargetVertexCurvature>& m_targetVertexCurvatures;
};
