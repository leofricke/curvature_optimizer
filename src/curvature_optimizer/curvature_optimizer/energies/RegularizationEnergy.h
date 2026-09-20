#pragma once

#include "opt/ModularGaussNewtonEnergyFunction.h"
#include <pmp/surface_mesh.h>

struct CurvatureOptimizerWorkspace;

struct RegularizationEnergy : public ModularGaussNewtonEnergyBlock
{
    struct InitialPosition
    {
        pmp::Vertex m_vertex;
        Eigen::Vector3d m_position;
    };

    RegularizationEnergy(const CurvatureOptimizerWorkspace& ws);

    void resetPositions();

    void evaluateResidual(const Eigen::VectorXd&, GaussNewtonSystem& system) override;
    void evaluateResidualAndJacobian(const Eigen::VectorXd&, GaussNewtonSystem& system) override;

    void initWeights(double energyWeight, Eigen::VectorXd& weights) override;

    const CurvatureOptimizerWorkspace& m_ws;
    std::vector<InitialPosition> m_initialPositions;
};
