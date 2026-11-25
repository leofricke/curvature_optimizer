#pragma once

#include "opt/LambdaEnergyFunction.h"

inline std::shared_ptr<GaussNewtonEnergyFunction> createRosenbrockFunctionResidual(size_t n)
{
    GaussNewtonEnergyLambdaFunction::ResidualLambda residualLambda = [n](const Eigen::VectorXd& x, GaussNewtonSystem& system)
    {
        for (size_t i = 0; i < n; i += 2)
        {
            system.m_residual[i] = 1.0 - x(i);
            system.m_residual[i + 1] = 10.0 * (x(i + 1) - x(i) * x(i));
        }
    };

    GaussNewtonEnergyLambdaFunction::ResidualAndJacobianLambda jacobianLambda = [n](const Eigen::VectorXd& x, GaussNewtonSystem& system)
    {
        system.beginFill();

        for (size_t i = 0; i < n; i += 2)
        {
            system.m_residual[i] = 1.0 - x(i);
            system.m_residual[i + 1] = 10.0 * (x(i + 1) - x(i) * x(i));
        }

        spmat::SpVec<double> row = system.createJacobianRow();

        for (size_t i = 0; i < n; i += 2)
        {
            row.clearPattern();
            row.beginFill();
            row.insert(i, -1.0);
            row.endFill();
            system.insertJacobianRow(i, row);

            row.clearPattern();
            row.beginFill();
            row.insert(i, -20.0 * x(i));
            row.insert(i + 1, 10.0);
            row.endFill();
            system.insertJacobianRow(i + 1, row);
        }

        system.endFill();
    };

    return std::make_shared<GaussNewtonEnergyLambdaFunction>(n, n, residualLambda, jacobianLambda);
}
