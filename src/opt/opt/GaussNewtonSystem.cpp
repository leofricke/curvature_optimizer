#include "GaussNewtonSystem.h"

GaussNewtonSystem::GaussNewtonSystem(size_t numResiduals, size_t numVars)
    : m_numResiduals(numResiduals), m_numVars(numVars), m_buildFullJacobian(false)
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
    m_hessianBuilder.setMatrix(&m_hessian);

    if (m_buildFullJacobian)
    {
        m_jacobian.resize(m_numResiduals, m_numVars);
    }
}

void GaussNewtonSystem::beginFilling()
{
    m_gradient.setZero();
    m_hessianBuilder.begin();

    if (m_buildFullJacobian)
    {
        m_jacobian.resize(m_numResiduals, m_numVars);
        m_jacobian.setZero();
    }
}

void GaussNewtonSystem::endFilling()
{
    m_hessianBuilder.end();
}

void GaussNewtonSystem::insertJacobianRow(Eigen::Index row, SparseVector& vec)
{
    const double r = m_residual[row];
    const double w = m_weights[row];

    vec.sortEntries();
    const std::vector<SparseVector::Entry>& entries = vec.getEntries();
    const size_t n = entries.size();

    for (size_t j = 0; j < n;)
    {
        const int col = entries[j].m_idx;

        double colValue = 0.0;
        size_t jEnd = j;
        while (jEnd < n && entries[jEnd].m_idx == col)
        {
            colValue += entries[jEnd].m_value;
            ++jEnd;
        }

        FixedPatternMatrixBuilder::ColumnInserter inserter = m_hessianBuilder.columnInserter(col);
        for (size_t i = j; i < n; ++i) // lower-triangular part only
        {
            const SparseVector::Entry& e0 = entries[i];
            inserter.addAscending(e0.m_idx, e0.m_value * w * colValue);
        }

        m_gradient[col] += colValue * w * r;

        if (m_buildFullJacobian)
        {
            m_jacobian.coeffRef(row, col) += colValue;
        }

        j = jEnd;
    }
}

double GaussNewtonSystem::evaluateEnergy() const
{
    return 0.5 * m_residual.dot(m_weights.asDiagonal() * m_residual);
}
