#pragma once

#include "GaussNewtonEnergyFunction.h"

/**
* Specialized version of GaussNewtonEnergyFunction that allows to provide residual and Jacobian evaluation via lambda functions.
*/
struct GaussNewtonEnergyLambdaFunction : public GaussNewtonEnergyFunction
{
    using ResidualLambda = std::function<void(const Eigen::VectorXd& params, GaussNewtonSystem& system)>;
    using ResidualAndJacobianLambda = std::function<void(const Eigen::VectorXd& params, GaussNewtonSystem& system)>;

    GaussNewtonEnergyLambdaFunction(size_t numResiduals, size_t numVars, ResidualLambda residualLambda, ResidualAndJacobianLambda residualAndJacobianLambda)
        : GaussNewtonEnergyFunction(numResiduals, numVars), m_residualLambda(residualLambda), m_residualAndJacobianLambda(residualAndJacobianLambda)
    {
    }

    ResidualLambda m_residualLambda;
    ResidualAndJacobianLambda m_residualAndJacobianLambda;

protected:
    double evaluateResidualImpl(const Eigen::VectorXd& params) override
    {
        m_residualLambda(params, m_system);
        return m_system.evaluateEnergy();
    }

    double evaluateResidualAndJacobianImpl(const Eigen::VectorXd& params) override
    {
        m_residualAndJacobianLambda(params, m_system);
        return m_system.evaluateEnergy();
    }
};