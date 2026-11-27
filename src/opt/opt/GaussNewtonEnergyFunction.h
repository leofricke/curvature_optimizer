#pragma once

#include "EnergyFunction.h"

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cassert>
#include <map>
#include <memory>
#include <spmat/SpMat.h>
#include <spmat/SpVec.h>
#include <string>
#include <vector>

struct GaussNewtonEnergyFunction;

using JacobianRow = spmat::SpVec<double>;

/**
* Access to the residuals and Jacobian rows of the element that is currently evaluated. Indices are local to the element.
*/
class GaussNewtonElementContext
{
public:
    bool withJacobian() const
    {
        return m_withJacobian;
    }

    const Eigen::VectorXd& params() const
    {
        return *m_params;
    }

    double& residual(size_t k);

    /**
    * Empty row for the derivatives of residual k. Rows are inserted after the element is evaluated.
    */
    JacobianRow& jacobianRow(size_t k);

    /**
    * Inserts the requested jacobian rows now.
    */
    void commitRows();

private:
    friend struct GaussNewtonEnergyFunction;

    GaussNewtonElementContext() = default;
    void beginEvaluation(GaussNewtonEnergyFunction& energyFunc, const Eigen::VectorXd& params, bool withJacobian);
    void beginElement(size_t offset, size_t numResiduals, const double* weights, double moduleWeight);

    GaussNewtonEnergyFunction* m_energyFunc = nullptr;
    const Eigen::VectorXd* m_params = nullptr;
    bool m_withJacobian = false;
    size_t m_offset = 0;
    size_t m_numResiduals = 0;
    const double* m_weights = nullptr;
    double m_moduleWeight = 1.0;
    std::vector<JacobianRow> m_rows;
    std::vector<bool> m_rowUsed;
};

/**
* One term of a GaussNewtonEnergyFunction: a weighted sum of squared residuals, grouped into elements (e.g. one per face).
*/
struct GaussNewtonEnergyModule
{
    GaussNewtonEnergyModule(std::string name, size_t numVars)
        : m_name(std::move(name)), m_numVars(numVars)
    {
    }

    virtual ~GaussNewtonEnergyModule() = default;

    virtual size_t numElements() const = 0;
    virtual size_t numResiduals(size_t element) const = 0;

    /**
    * Computes the residuals (and Jacobian rows if ctx.withJacobian()) of one element.
    * Called concurrently for different elements: must only write to ctx.
    */
    virtual void evaluateElement(size_t element, GaussNewtonElementContext& ctx) const = 0;

    /**
    * Called once the residual layout is known. Modules can set their weights here.
    */
    virtual void initialize()
    {
    }

    /**
    * Weights of all residuals of the module, in element order. Must not be changed during an evaluation.
    */
    Eigen::Ref<Eigen::VectorXd> residualWeights()
    {
        return m_residualWeights;
    }

    const Eigen::VectorXd& residualWeights() const
    {
        return m_residualWeights;
    }

    /**
    * Weights of the residuals of one element. Available from initialize() on.
    */
    Eigen::VectorBlock<Eigen::VectorXd> residualWeights(size_t element)
    {
        assert(element + 1 < m_elementOffsets.size() && "residualWeights(element) is only available after initialize() was called");
        return m_residualWeights.segment(static_cast<Eigen::Index>(m_elementOffsets[element]), static_cast<Eigen::Index>(m_elementOffsets[element + 1] - m_elementOffsets[element]));
    }

    std::string m_name;
    size_t m_numVars;

private:
    friend struct GaussNewtonEnergyFunction;

    Eigen::VectorXd m_residualWeights;
    std::vector<size_t> m_elementOffsets; // residual offset of each element within the module, size numElements() + 1
    size_t m_offset = 0;                  // first residual of the module in the energy function
    double m_moduleWeight = 1.0;
    const GaussNewtonEnergyFunction* m_owner = nullptr; // a module can only be part of one energy function
};

/**
* Energy E(x) = 0.5 * sum_i w_i * r_i(x)^2, given as the sum of weighted modules.
* Builds the Gauss-Newton system J^T W J * delta = -J^T W r, evaluating all elements of all modules in parallel.
*/
struct GaussNewtonEnergyFunction : public EnergyFunction
{
    using HessianType = spmat::SpMat<double, spmat::SpMatLayout::ColMajor, spmat::SpMatSymmetry::Lower, spmat::SpMatThreadSafety::Enabled>;

    /**
    * moduleWeights scales the energy of each module (all ones if empty).
    */
    GaussNewtonEnergyFunction(const std::vector<std::shared_ptr<GaussNewtonEnergyModule>>& modules, const std::vector<double>& moduleWeights = {});
    ~GaussNewtonEnergyFunction() override;

    GaussNewtonEnergyFunction(const GaussNewtonEnergyFunction&) = delete;
    GaussNewtonEnergyFunction& operator=(const GaussNewtonEnergyFunction&) = delete;

    const Eigen::VectorXd& gradient() const override
    {
        return m_gradient;
    }

    /**
    * Weighted energies from the residuals of the last evaluation.
    */
    double getEnergy() const;
    double getModuleEnergy(size_t idx) const;
    std::map<std::string, double> getModuleEnergies() const;

    std::vector<std::shared_ptr<GaussNewtonEnergyModule>> m_modules;

    size_t m_numResiduals;
    Eigen::VectorXd m_residual;
    Eigen::VectorXd m_gradient; // J^T * W * r
    HessianType m_hessian;      // J^T * W * J
    bool m_hessianPatternDirty; // set if the last evaluation changed the sparsity pattern of the Hessian

    bool m_buildFullJacobian; // full jacobian only for unit tests, too expensive otherwise
    spmat::SpMat<double, spmat::SpMatLayout::RowMajor, spmat::SpMatSymmetry::None, spmat::SpMatThreadSafety::Enabled> m_jacobian;

protected:
    double evaluateEnergyImpl(const Eigen::VectorXd& params) override;
    double evaluateEnergyAndGradientImpl(const Eigen::VectorXd& params) override;

private:
    friend class GaussNewtonElementContext;

    void evaluateElements(const Eigen::VectorXd& params, bool withJacobian);
    void insertJacobianRow(size_t row, const JacobianRow& vec, double weight);

    std::vector<GaussNewtonElementContext> m_contexts;     // one per thread, kept to reuse the row allocations
    std::vector<size_t> m_moduleElementOffsets; // first element of each module among all elements, size numModules + 1
};
