#pragma once

#include "opt/LambdaEnergyFunction.h"

// (x - 2)^2 + 2(y + 3)^2 + 3(z - 5)^2
inline std::shared_ptr<GaussNewtonEnergyFunction> createSimpleQuadraticFunctionResidual()
{
    GaussNewtonEnergyLambdaFunction::ResidualLambda residualLambda = [](const Eigen::VectorXd& x, GaussNewtonSystem& system)
    {
        system.m_residual[0] = std::sqrt(2.0) * (x[0] - 2.0);
        system.m_residual[1] = std::sqrt(4.0) * (x[1] + 3.0);
        system.m_residual[2] = std::sqrt(6.0) * (x[2] - 5.0);
    };

    GaussNewtonEnergyLambdaFunction::ResidualAndJacobianLambda jacobianLambda = [](const Eigen::VectorXd& x, GaussNewtonSystem& system)
    {
        system.beginFill();

        system.m_residual[0] = std::sqrt(2.0) * (x[0] - 2.0);
        system.m_residual[1] = std::sqrt(4.0) * (x[1] + 3.0);
        system.m_residual[2] = std::sqrt(6.0) * (x[2] - 5.0);

        spmat::SpVec<double> row = system.createJacobianRow();

        row.clearPattern();
        row.beginFill();
        row.insert(0, std::sqrt(2.0));
        row.endFill();
        system.insertJacobianRow(0, row);

        row.clearPattern();
        row.beginFill();
        row.insert(1, std::sqrt(4.0));
        row.endFill();
        system.insertJacobianRow(1, row);

        row.clearPattern();
        row.beginFill();
        row.insert(2, std::sqrt(6.0));
        row.endFill();
        system.insertJacobianRow(2, row);

        system.endFill();
    };

    return std::make_shared<GaussNewtonEnergyLambdaFunction>(3, 3, residualLambda, jacobianLambda);
}
