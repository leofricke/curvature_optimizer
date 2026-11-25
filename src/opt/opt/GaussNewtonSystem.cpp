#include "GaussNewtonSystem.h"

GaussNewtonSystem::GaussNewtonSystem(size_t numResiduals, size_t numVars)
    : m_numResiduals(numResiduals), m_numVars(numVars), m_hessianPatternDirty(true), m_buildFullJacobian(false)
{
    init(numResiduals, numVars);
}

void GaussNewtonSystem::init(size_t numResiduals, size_t numVars)
{
    m_numResiduals = numResiduals;
    m_numVars = numVars;

    m_residual.resize(numResiduals);
    m_weights.resize(numResiduals);
    m_weights.setOnes();
    m_gradient.resize(numVars);
    m_hessian.resize(numVars, numVars);
    m_hessian.clearPattern();
    m_hessianPatternDirty = true;

    if (m_buildFullJacobian)
    {
        m_jacobian.resize(m_numResiduals, m_numVars);
    }
}

void GaussNewtonSystem::beginFill()
{
    m_gradient.setZero();
    m_hessian.beginFill();
    m_hessianPatternDirty = false;

    if (m_hessian.nonZeros() == 0)
    {
        // create explicit diagonal
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(m_numVars); ++i)
        {
            m_hessian.insert(i, i, 0.0);
        }
    }

    if (m_buildFullJacobian)
    {
        m_jacobian.resize(m_numResiduals, m_numVars);
        m_jacobian.setZero();
    }
}

bool GaussNewtonSystem::endFill()
{
    m_hessianPatternDirty = m_hessian.endFill();
    return m_hessianPatternDirty;
}

void GaussNewtonSystem::insertJacobianRow(Eigen::Index row, const spmat::SpVec<double>& vec)
{
    const double r = m_residual[row];
    const double w = m_weights[row];

    for (const spmat::SpVec<double>::Entry& e : vec.entries())
    {
        m_gradient[e.m_idx] += e.m_value * w * r;
    }
    m_hessian.insertOuterProduct(vec, w);

    if (m_buildFullJacobian)
    {
        for (const spmat::SpVec<double>::Entry& e : vec.entries())
        {
            m_jacobian.coeffRef(row, e.m_idx) += e.m_value;
        }
    }
}

spmat::SpVec<double> GaussNewtonSystem::createJacobianRow() const
{
    return spmat::SpVec<double>(static_cast<int>(m_numVars));
}

double GaussNewtonSystem::evaluateEnergy() const
{
    return 0.5 * m_residual.dot(m_weights.asDiagonal() * m_residual);
}
