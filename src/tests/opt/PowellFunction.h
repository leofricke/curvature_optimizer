#pragma once

#include "opt/LambdaEnergyFunction.h"

inline std::shared_ptr<GaussNewtonEnergyFunction> createPowellFunctionResidual()
{
    GaussNewtonEnergyLambdaFunction::ResidualLambda residualLambda = [](const Eigen::VectorXd& x, GaussNewtonSystem& system)
    {
        system.m_residual[0] = x[0] + 10.0 * x[1];
        system.m_residual[1] = std::sqrt(5.0) * (x[2] - x[3]);
        system.m_residual[2] = std::pow(x[1] - 2.0 * x[2], 2.0);
        system.m_residual[3] = std::sqrt(10.0) * std::pow(x[0] - x[3], 2.0);
    };

    GaussNewtonEnergyLambdaFunction::ResidualAndJacobianLambda jacobianLambda = [](const Eigen::VectorXd& x, GaussNewtonSystem& system)
    {
        system.beginFill();

        system.m_residual[0] = x[0] + 10.0 * x[1];
        system.m_residual[1] = std::sqrt(5.0) * (x[2] - x[3]);
        system.m_residual[2] = std::pow(x[1] - 2.0 * x[2], 2.0);
        system.m_residual[3] = std::sqrt(10.0) * std::pow(x[0] - x[3], 2.0);

        spmat::SpVec<double> row = system.createJacobianRow();

        row.clearPattern();
        row.beginFill();
        row.insert(0, 1.0);
        row.insert(1, 10.0);
        row.endFill();
        system.insertJacobianRow(0, row);

        row.clearPattern();
        row.beginFill();
        row.insert(2, std::sqrt(5.0));
        row.insert(3, -std::sqrt(5.0));
        row.endFill();
        system.insertJacobianRow(1, row);

        row.clearPattern();
        row.beginFill();
        row.insert(1, 2.0 * (x[1] - 2.0 * x[2]));
        row.insert(2, -4.0 * (x[1] - 2.0 * x[2]));
        row.endFill();
        system.insertJacobianRow(2, row);

        row.clearPattern();
        row.beginFill();
        row.insert(0, 2.0 * std::sqrt(10.0) * (x[0] - x[3]));
        row.insert(3, -2.0 * std::sqrt(10.0) * (x[0] - x[3]));
        row.endFill();
        system.insertJacobianRow(3, row);

        system.endFill();
    };

    return std::make_shared<GaussNewtonEnergyLambdaFunction>(4, 4, residualLambda, jacobianLambda);
}
