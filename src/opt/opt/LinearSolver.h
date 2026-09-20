#pragma once

#include <memory>
#include <Eigen/Sparse>

class LinearSolverImpl
{
public:
    LinearSolverImpl() = default;
    virtual ~LinearSolverImpl() = default;

    virtual bool compute(const Eigen::SparseMatrix<double, Eigen::ColMajor>& A, bool patternDirty) = 0;
    virtual bool solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution) = 0;
};

class LinearSolver
{
public:
    LinearSolver();

    bool compute(const Eigen::SparseMatrix<double, Eigen::ColMajor>& A, bool patternDirty);
    bool solve(const Eigen::VectorXd& rhs, Eigen::VectorXd& solution);

private:
    std::unique_ptr<LinearSolverImpl> m_impl;
};
