#pragma once

#include "utils/SparseVector.h"
#include "utils/FixedPatternMatrixBuilder.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

/*
* Represents the Gauss-Newton system J(x)^T * W * J(x) * delta = -J(x)^T * W * r(x),
* where J(x) is the Jacobian of r(x) and W is a diagonal matrix containing residual weights.
* J(x)^T * W * r(x) is the gradient of the energy function, and J(x)^T * W * J(x) is the Hessian replacement.
* To make optimization efficient, this class builds the system without explicitly forming the Jacobian matrix,
* but instead directly accumulates the Hessian replacement and gradient from individual Jacobian rows.
* The full Jacobian may be optionally requested for testing purposes.
*/
struct GaussNewtonSystem
{
    /**
    * Initialize a system of given size.
    */
    GaussNewtonSystem(size_t numResiduals, size_t numVars);

    /**
    * Initialize a system of given size.
    */
    void init(size_t numResiduals, size_t numVars);
    /**
    * Must be called before inserting Jacobian rows. Resets the gradient and Hessian to zero.
    */
    void beginFilling();
    /**
    * Must be called after all Jacobian rows are inserted.
    */
    void endFilling();

    /**
    * Inserts one row of the Jacobian into the system.
    */
    void insertJacobianRow(Eigen::Index row, SparseVector& vec);
    /**
    * Helper to evaluate 0.5 * r^T * W * r.
    */
    double evaluateEnergy() const;

    size_t m_numResiduals;
    size_t m_numVars;
    Eigen::VectorXd m_residual;
    Eigen::VectorXd m_weights;
    Eigen::VectorXd m_gradient;            // J^T * W * r
    Eigen::SparseMatrix<double> m_hessian; // J^T * W * J
    FixedPatternMatrixBuilder m_hessianBuilder;

    bool m_buildFullJacobian; // use only for unit tests, otherwise this is too expensive
    Eigen::SparseMatrix<double> m_jacobian;
};
