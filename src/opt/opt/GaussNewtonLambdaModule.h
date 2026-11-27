#pragma once

#include "GaussNewtonEnergyFunction.h"

#include <functional>

/**
* Energy module defined by lambda functions.
*/
struct GaussNewtonLambdaModule : public GaussNewtonEnergyModule
{
    using NumResidualsLambda = std::function<size_t(size_t element)>;
    using EvaluateElementLambda = std::function<void(size_t element, GaussNewtonElementContext& ctx)>;

    GaussNewtonLambdaModule(std::string name, size_t numVars, size_t numElements, NumResidualsLambda numResidualsLambda, EvaluateElementLambda evaluateElementLambda)
        : GaussNewtonEnergyModule(std::move(name), numVars), m_numElements(numElements), m_numResidualsLambda(std::move(numResidualsLambda)), m_evaluateElementLambda(std::move(evaluateElementLambda))
    {
    }

    size_t numElements() const override
    {
        return m_numElements;
    }

    size_t numResiduals(size_t element) const override
    {
        return m_numResidualsLambda(element);
    }

    void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const override
    {
        m_evaluateElementLambda(element, ctx);
    }

    size_t m_numElements;
    NumResidualsLambda m_numResidualsLambda;
    EvaluateElementLambda m_evaluateElementLambda;
};
