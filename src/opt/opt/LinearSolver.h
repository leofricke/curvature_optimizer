#pragma once

#include <memory>
#include <Eigen/Sparse>

/**
* Read-only view of a column-major sparse matrix.
* Only the lower triangular part of the symmetric system matrix is used.
*/
using SparseMatrixView = Eigen::Map<const Eigen::SparseMatrix<double, Eigen::ColMajor, int>>;

class LinearSolverImpl
{
public:
    LinearSolverImpl() = default;
    virtual ~LinearSolverImpl() = default;

    virtual bool compute(const SparseMatrixView& A, bool patternDirty) = 0;
    virtual bool solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution) = 0;
};

class LinearSolver
{
public:
    LinearSolver();

    bool compute(const SparseMatrixView& A, bool patternDirty);
    bool solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution);

private:
    std::unique_ptr<LinearSolverImpl> m_impl;
};
