#pragma once

#include "opt/GaussNewtonLambdaModule.h"

inline std::shared_ptr<GaussNewtonEnergyFunction> createPowellFunctionResidual()
{
    auto module = std::make_shared<GaussNewtonLambdaModule>(
        "Powell", 4, 1, [](size_t /*element*/) { return size_t(4); },
        [](size_t /*element*/, GaussNewtonElementContext& ctx)
        {
            const Eigen::VectorXd& x = ctx.params();
            ctx.residual(0) = x[0] + 10.0 * x[1];
            ctx.residual(1) = std::sqrt(5.0) * (x[2] - x[3]);
            ctx.residual(2) = std::pow(x[1] - 2.0 * x[2], 2.0);
            ctx.residual(3) = std::sqrt(10.0) * std::pow(x[0] - x[3], 2.0);
            if (ctx.withJacobian())
            {
                ctx.jacobianRow(0).insert(0, 1.0);
                ctx.jacobianRow(0).insert(1, 10.0);
                ctx.jacobianRow(1).insert(2, std::sqrt(5.0));
                ctx.jacobianRow(1).insert(3, -std::sqrt(5.0));
                ctx.jacobianRow(2).insert(1, 2.0 * (x[1] - 2.0 * x[2]));
                ctx.jacobianRow(2).insert(2, -4.0 * (x[1] - 2.0 * x[2]));
                ctx.jacobianRow(3).insert(0, 2.0 * std::sqrt(10.0) * (x[0] - x[3]));
                ctx.jacobianRow(3).insert(3, -2.0 * std::sqrt(10.0) * (x[0] - x[3]));
            }
        });
    return std::make_shared<GaussNewtonEnergyFunction>(std::vector<std::shared_ptr<GaussNewtonEnergyModule>> { module });
}
