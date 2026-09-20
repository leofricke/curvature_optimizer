#pragma once

#include "EnergyFunction.h"
#include "GaussNewtonSystem.h"
#include <Eigen/Core>

/**
* Gauss-Newton energy function.
* Inherits from EnergyFunction and models an energy in the form E(x) = 0.5 * \sum_i ||r_i(x)||^2, where r_i(x) are the residuals.
* Residuals are stored in the residual vector r(x) such that E(x) = 0.5 * r(x)^T * r(x).
* For Gauss-Newton optimization, we must build the system J(x)^T * W * J(x) * delta = -J(x)^T * W * r(x), where J(x) is the Jacobian of r(x).
*/
struct GaussNewtonEnergyFunction : public EnergyFunction
{
    GaussNewtonEnergyFunction(size_t numResiduals, size_t numVars)
        : EnergyFunction(numVars), m_numResiduals(numResiduals), m_system(numResiduals, numVars)
    {
    }

    virtual ~GaussNewtonEnergyFunction() = default;

    const Eigen::VectorXd& gradient() const override
    {
        return m_system.m_gradient;
    }

    double evaluateResidual(const Eigen::VectorXd& params)
    {
        invokeBeforeEvaluateCallback(params, DerivativeMode::NoDerivs);
        return evaluateResidualImpl(params);
    }

    double evaluateResidualAndJacobian(const Eigen::VectorXd& params)
    {
        invokeBeforeEvaluateCallback(params, DerivativeMode::FirstDerivs);
        return evaluateResidualAndJacobianImpl(params);
    }

    size_t m_numResiduals;
    GaussNewtonSystem m_system;

protected:
    virtual double evaluateResidualImpl(const Eigen::VectorXd& params) = 0;
    virtual double evaluateResidualAndJacobianImpl(const Eigen::VectorXd& params) = 0;

    double evaluateEnergyImpl(const Eigen::VectorXd& params) override
    {
        return evaluateResidualImpl(params);
    }

    double evaluateEnergyAndGradientImpl(const Eigen::VectorXd& params) override
    {
        return evaluateResidualAndJacobianImpl(params);
    }
};
