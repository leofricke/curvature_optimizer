#pragma once

#include "opt/GaussNewtonLambdaModule.h"

inline std::shared_ptr<GaussNewtonEnergyFunction> createRosenbrockFunctionResidual(size_t n)
{
    auto module = std::make_shared<GaussNewtonLambdaModule>(
        "Rosenbrock", n, n / 2, [](size_t /*element*/) { return size_t(2); },
        [](size_t element, GaussNewtonElementContext& ctx)
        {
            const Eigen::VectorXd& x = ctx.params();
            const int i = static_cast<int>(element) * 2;
            ctx.residual(0) = 1.0 - x(i);
            ctx.residual(1) = 10.0 * (x(i + 1) - x(i) * x(i));
            if (ctx.withJacobian())
            {
                ctx.jacobianRow(0).insert(i, -1.0);
                ctx.jacobianRow(1).insert(i, -20.0 * x(i));
                ctx.jacobianRow(1).insert(i + 1, 10.0);
            }
        });
    return std::make_shared<GaussNewtonEnergyFunction>(std::vector<std::shared_ptr<GaussNewtonEnergyModule>> { module });
}
