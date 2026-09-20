#pragma once

#include <Eigen/Core>
#include <functional>
#include <cassert>
#include <vector>

/**
* Base class for an energy function.
* Provides an interface for evaluating the energy and its gradient.
* Also supports registering callbacks that run before evaluation.
*/
struct EnergyFunction
{
    enum class DerivativeMode
    {
        NoDerivs,
        FirstDerivs
        //SecondDerivs
    };
    using BeforeEvaluateCallback = std::function<void(const Eigen::VectorXd& params, DerivativeMode mode)>;

    EnergyFunction(size_t numVars)
        : m_numVars(numVars)
    {
    }

    virtual ~EnergyFunction() = default;

    virtual const Eigen::VectorXd& gradient() const = 0;

    /**
    * Evaluates the energy function at parameter vector params.
    */
    double evaluateEnergy(const Eigen::VectorXd& params)
    {
        assert(static_cast<size_t>(params.size()) == m_numVars);
        invokeBeforeEvaluateCallback(params, DerivativeMode::NoDerivs);
        return evaluateEnergyImpl(params);
    }

    /**
    * Evaluates the energy function and its gradient at parameter vector params.
    * After calling this function, the gradient can be accessed via the gradient() method.
    */
    double evaluateEnergyAndGradient(const Eigen::VectorXd& params)
    {
        assert(static_cast<size_t>(params.size()) == m_numVars);
        invokeBeforeEvaluateCallback(params, DerivativeMode::FirstDerivs);
        return evaluateEnergyAndGradientImpl(params);
    }

    void addBeforeEvaluateCallback(const BeforeEvaluateCallback& callback)
    {
        m_beforeEvaluateCallbacks.push_back(callback);
    }

    size_t m_numVars;

protected:
    std::vector<BeforeEvaluateCallback> m_beforeEvaluateCallbacks;

    void invokeBeforeEvaluateCallback(const Eigen::VectorXd& params, DerivativeMode mode) const
    {
        for (const auto& callback : m_beforeEvaluateCallbacks)
        {
            callback(params, mode);
        }
    }

    virtual double evaluateEnergyImpl(const Eigen::VectorXd& params) = 0;
    virtual double evaluateEnergyAndGradientImpl(const Eigen::VectorXd& params) = 0;
};