#pragma once

#include "opt/GaussNewtonLambdaModule.h"

// (x - 2)^2 + 2(y + 3)^2 + 3(z - 5)^2
inline std::shared_ptr<GaussNewtonEnergyFunction> createSimpleQuadraticFunctionResidual()
{
    auto module = std::make_shared<GaussNewtonLambdaModule>(
        "SimpleQuadratic", 3, 3, [](size_t /*element*/) { return size_t(1); },
        [](size_t element, GaussNewtonElementContext& ctx)
        {
            const Eigen::VectorXd& x = ctx.params();
            const int i = static_cast<int>(element);
            const double scale = std::sqrt(2.0 * (i + 1));
            const double target[3] = { 2.0, -3.0, 5.0 };
            ctx.residual(0) = scale * (x[i] - target[i]);
            if (ctx.withJacobian())
            {
                ctx.jacobianRow(0).insert(i, scale);
            }
        });
    return std::make_shared<GaussNewtonEnergyFunction>(std::vector<std::shared_ptr<GaussNewtonEnergyModule>> { module });
}
