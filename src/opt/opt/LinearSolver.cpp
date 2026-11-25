#include "LinearSolver.h"

class EigenLDLTSolver : public LinearSolverImpl
{
public:
    EigenLDLTSolver() = default;

    bool compute(const SparseMatrixView& A, bool patternDirty) override
    {
        if (patternDirty)
        {
            m_solver.analyzePattern(A);
        }
        m_solver.factorize(A);
        return m_solver.info() == Eigen::Success;
    }

    bool solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution) override
    {
        solution = m_solver.solve(rhs);
        return m_solver.info() == Eigen::Success;
    }

private:
    Eigen::SimplicialLDLT<SparseMatrixView, Eigen::Lower> m_solver;
};

#if defined(HAS_MKL)
#include <Eigen/PardisoSupport>
#endif
class PardisoLDLTSolver : public LinearSolverImpl
{
public:
    PardisoLDLTSolver()
    {
        //m_solver.pardisoParameterArray()[33] = omp_get_num_threads(); // enable CNR mode for reproducible results
    }

    bool compute([[maybe_unused]] const SparseMatrixView& A, [[maybe_unused]] bool patternDirty) override
    {
#if defined(HAS_MKL)
        if (patternDirty)
        {
            m_solver.analyzePattern(A);
        }
        m_solver.factorize(A);
        return m_solver.info() == Eigen::Success;
#else
        return false;
#endif
    }

    bool solve([[maybe_unused]] const Eigen::VectorXd& rhs, [[maybe_unused]] Eigen::VectorXd& solution) override
    {
#if defined(HAS_MKL)
        solution = m_solver.solve(rhs);
        return m_solver.info() == Eigen::Success;
#else
        return false;
#endif
    }

private:
#if defined(HAS_MKL)
    Eigen::PardisoLDLT<SparseMatrixView, Eigen::Lower> m_solver;
#endif
};

#if defined(HAS_ACCELERATE)
#include <Eigen/AccelerateSupport>
#endif
class AccelerateLDLTSolver : public LinearSolverImpl
{
public:
    AccelerateLDLTSolver() = default;

    bool compute([[maybe_unused]] const SparseMatrixView& A, [[maybe_unused]] bool patternDirty) override
    {
#if defined(HAS_ACCELERATE)
        if (patternDirty)
        {
            m_solver.analyzePattern(A);
        }
        m_solver.factorize(A);
        return m_solver.info() == Eigen::Success;
#else
        return false;
#endif
    }

    bool solve([[maybe_unused]] const Eigen::VectorXd& rhs, [[maybe_unused]] Eigen::VectorXd& solution) override
    {
#if defined(HAS_ACCELERATE)
        solution = m_solver.solve(rhs);
        return m_solver.info() == Eigen::Success;
#else
        return false;
#endif
    }

private:
#if defined(HAS_ACCELERATE)
    Eigen::AccelerateLDLT<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::Lower> m_solver;
#endif
};

LinearSolver::LinearSolver()
{
#if defined(HAS_MKL)
    m_impl = std::make_unique<PardisoLDLTSolver>();
#elif defined(HAS_ACCELERATE)
    m_impl = std::make_unique<AccelerateLDLTSolver>();
#else
    m_impl = std::make_unique<EigenLDLTSolver>();
#endif
}

bool LinearSolver::compute(const SparseMatrixView& A, bool patternDirty)
{
    return m_impl->compute(A, patternDirty);
}

bool LinearSolver::solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution)
{
    return m_impl->solve(rhs, solution);
}
