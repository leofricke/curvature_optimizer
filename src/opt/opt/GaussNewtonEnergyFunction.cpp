#include "GaussNewtonEnergyFunction.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <set>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

double& GaussNewtonElementContext::residual(size_t k)
{
    assert(k < m_numResiduals);
    return m_energyFunc->m_residual[m_offset + k];
}

JacobianRow& GaussNewtonElementContext::jacobianRow(size_t k)
{
    assert(m_withJacobian && "jacobianRow() requires withJacobian()");
    assert(k < m_numResiduals);
    if (!m_rowUsed[k])
    {
        m_rows[k].clearPattern();
        m_rows[k].beginFill();
        m_rowUsed[k] = true;
    }
    return m_rows[k];
}

void GaussNewtonElementContext::commitRows()
{
    const size_t n = std::min(m_numResiduals, m_rows.size());
    for (size_t k = 0; k < n; ++k)
    {
        if (m_rowUsed[k])
        {
            m_rows[k].endFill();
            m_energyFunc->insertJacobianRow(m_offset + k, m_rows[k], m_moduleWeight * m_weights[k]);
            m_rowUsed[k] = false;
        }
    }
}

void GaussNewtonElementContext::beginEvaluation(GaussNewtonEnergyFunction& energyFunc, const Eigen::VectorXd& params, bool withJacobian)
{
    m_energyFunc = &energyFunc;
    m_params = &params;
    m_withJacobian = withJacobian;
}

void GaussNewtonElementContext::beginElement(size_t offset, size_t numResiduals, const double* weights, double moduleWeight)
{
    m_offset = offset;
    m_numResiduals = numResiduals;
    if (m_withJacobian && numResiduals > m_rows.size())
    {
        m_rows.resize(numResiduals, JacobianRow(static_cast<int>(m_energyFunc->m_numVars)));
        m_rowUsed.resize(numResiduals, false);
    }
    m_weights = weights;
    m_moduleWeight = moduleWeight;
}

GaussNewtonEnergyFunction::GaussNewtonEnergyFunction(const std::vector<std::shared_ptr<GaussNewtonEnergyModule>>& modules, const std::vector<double>& moduleWeights)
    : EnergyFunction(0), m_modules(modules), m_numResiduals(0), m_hessianPatternDirty(true), m_buildFullJacobian(false)
{
    if (m_modules.empty())
    {
        throw std::runtime_error("GaussNewtonEnergyFunction requires at least one module");
    }
    if (!moduleWeights.empty() && moduleWeights.size() != m_modules.size())
    {
        throw std::runtime_error("Weights array has wrong size");
    }

    std::set<std::string> names;
    for (const auto& module : m_modules)
    {
        if (module->m_owner)
        {
            throw std::runtime_error("Module '" + module->m_name + "' is already part of another energy function");
        }
        if (!names.insert(module->m_name).second)
        {
            throw std::runtime_error("Module names must be unique: '" + module->m_name + "'");
        }
        if (module->m_numVars != m_modules[0]->m_numVars)
        {
            throw std::runtime_error("Number of variables must be equal for all modules");
        }
    }

    m_numVars = m_modules[0]->m_numVars;
    m_moduleElementOffsets.assign(1, 0);
    for (size_t m = 0; m < m_modules.size(); ++m)
    {
        GaussNewtonEnergyModule& module = *m_modules[m];
        module.m_owner = this;

        const size_t numElements = module.numElements();
        module.m_offset = m_numResiduals;
        module.m_moduleWeight = moduleWeights.empty() ? 1.0 : moduleWeights[m];
        module.m_elementOffsets.resize(numElements + 1);
        module.m_elementOffsets[0] = 0;
        for (size_t e = 0; e < numElements; ++e)
        {
            module.m_elementOffsets[e + 1] = module.m_elementOffsets[e] + module.numResiduals(e);
        }

        module.m_residualWeights.setOnes(static_cast<Eigen::Index>(module.m_elementOffsets[numElements]));
        module.initialize();

        m_numResiduals += module.m_elementOffsets[numElements];
        m_moduleElementOffsets.push_back(m_moduleElementOffsets.back() + numElements);
    }

    m_residual.resize(static_cast<Eigen::Index>(m_numResiduals));
    m_gradient.resize(static_cast<Eigen::Index>(m_numVars));
    m_hessian.resize(static_cast<int>(m_numVars), static_cast<int>(m_numVars));
}

GaussNewtonEnergyFunction::~GaussNewtonEnergyFunction()
{
    for (const auto& module : m_modules)
    {
        module->m_owner = nullptr;
    }
}

double GaussNewtonEnergyFunction::getEnergy() const
{
    double energy = 0.0;
    for (size_t i = 0; i < m_modules.size(); ++i)
    {
        energy += getModuleEnergy(i);
    }
    return energy;
}

double GaussNewtonEnergyFunction::getModuleEnergy(size_t idx) const
{
    const GaussNewtonEnergyModule& module = *m_modules[idx];
    const auto r = m_residual.segment(static_cast<Eigen::Index>(module.m_offset), module.m_residualWeights.size());
    return 0.5 * module.m_moduleWeight * r.dot(module.m_residualWeights.asDiagonal() * r);
}

std::map<std::string, double> GaussNewtonEnergyFunction::getModuleEnergies() const
{
    std::map<std::string, double> energies;
    for (size_t i = 0; i < m_modules.size(); ++i)
    {
        energies.emplace(m_modules[i]->m_name, getModuleEnergy(i));
    }
    return energies;
}

double GaussNewtonEnergyFunction::evaluateEnergyImpl(const Eigen::VectorXd& params)
{
    evaluateElements(params, false);
    return getEnergy();
}

double GaussNewtonEnergyFunction::evaluateEnergyAndGradientImpl(const Eigen::VectorXd& params)
{
    m_gradient.setZero();
    m_hessian.beginFill();
    if (m_hessian.nonZeros() == 0)
    {
        // explicit diagonal, needed for regularization
        for (int i = 0; i < static_cast<int>(m_numVars); ++i)
        {
            m_hessian.insert(i, i, 0.0);
        }
    }
    if (m_buildFullJacobian)
    {
        m_jacobian.resize(static_cast<int>(m_numResiduals), static_cast<int>(m_numVars));
        m_jacobian.beginFill();
    }

    evaluateElements(params, true);

    m_hessianPatternDirty = m_hessian.endFill();
    if (m_buildFullJacobian)
    {
        m_jacobian.endFill();
    }
    return getEnergy();
}

void GaussNewtonEnergyFunction::evaluateElements(const Eigen::VectorXd& params, bool withJacobian)
{
    const int numElements = static_cast<int>(m_moduleElementOffsets.back());
#ifdef _OPENMP
    const size_t numThreads = static_cast<size_t>(omp_get_max_threads());
#else
    const size_t numThreads = 1;
#endif
    if (m_contexts.size() < numThreads)
    {
        m_contexts.resize(numThreads, GaussNewtonElementContext());
    }

#pragma omp parallel
    {
#ifdef _OPENMP
        GaussNewtonElementContext& ctx = m_contexts[static_cast<size_t>(omp_get_thread_num())];
#else
        GaussNewtonElementContext& ctx = m_contexts[0];
#endif
        ctx.beginEvaluation(*this, params, withJacobian);
#pragma omp for schedule(dynamic, 64)
        for (int element = 0; element < numElements; ++element)
        {
            const size_t m = static_cast<size_t>(std::upper_bound(m_moduleElementOffsets.begin(), m_moduleElementOffsets.end(), static_cast<size_t>(element)) - m_moduleElementOffsets.begin()) - 1;
            const GaussNewtonEnergyModule& module = *m_modules[m];
            const size_t localElement = static_cast<size_t>(element) - m_moduleElementOffsets[m];
            const size_t localOffset = module.m_elementOffsets[localElement];

            ctx.beginElement(module.m_offset + localOffset, module.m_elementOffsets[localElement + 1] - localOffset, module.m_residualWeights.data() + localOffset, module.m_moduleWeight);
            module.evaluateElement(localElement, ctx);
            if (withJacobian)
            {
                ctx.commitRows();
            }
        }
    }
}

void GaussNewtonEnergyFunction::insertJacobianRow(size_t row, const JacobianRow& vec, double weight)
{
    const double r = m_residual[static_cast<Eigen::Index>(row)];
    for (const JacobianRow::Entry& e : vec.entries())
    {
        std::atomic_ref<double>(m_gradient[e.m_idx]).fetch_add(e.m_value * weight * r, std::memory_order_relaxed);
    }
    m_hessian.insertOuterProduct(vec, weight);

    if (m_buildFullJacobian)
    {
        for (const JacobianRow::Entry& e : vec.entries())
        {
            m_jacobian.insert(static_cast<int>(row), e.m_idx, e.m_value);
        }
    }
}
