#include <gtest/gtest.h>

#include "PowellFunction.h"
#include "RosenbrockFunction.h"
#include "SimpleQuadraticFunction.h"
#include "opt/GaussNewtonOptimizer.h"

TEST(GaussNewtonOptimizerTest, Powell)
{
    std::shared_ptr<GaussNewtonEnergyFunction> energyFunc = createPowellFunctionResidual();
    GaussNewtonOptimizer optimizer;
    optimizer.setEnergyReductionTolerance(0.0);

    Eigen::VectorXd x = Eigen::VectorXd::Constant(4, 10.0);
    const double finalEnergy = optimizer.run(energyFunc, x);

    for (size_t i = 0; i < 4; ++i)
    {
        EXPECT_NEAR(x[i], 0.0, 1e-3);
    }
    EXPECT_NEAR(finalEnergy, 0.0, 1e-8);
}

TEST(GaussNewtonOptimizerTest, Rosenbrock100)
{
    const size_t n = 100;
    std::shared_ptr<GaussNewtonEnergyFunction> energyFunc = createRosenbrockFunctionResidual(n);
    GaussNewtonOptimizer optimizer;
    optimizer.setEnergyReductionTolerance(0.0);

    Eigen::VectorXd x = Eigen::VectorXd::Constant(n, 3.0);
    const double finalEnergy = optimizer.run(energyFunc, x);

    for (size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(x[i], 1.0, 1e-8);
    }
    EXPECT_NEAR(finalEnergy, 0.0, 1e-8);
}

TEST(GaussNewtonOptimizerTest, SimpleQuadratic)
{
    std::shared_ptr<GaussNewtonEnergyFunction> energyFunc = createSimpleQuadraticFunctionResidual();
    GaussNewtonOptimizer optimizer;
    optimizer.setEnergyReductionTolerance(0.0);

    Eigen::VectorXd x = Eigen::VectorXd::Constant(3, 42.0);
    const double finalEnergy = optimizer.run(energyFunc, x);

    EXPECT_NEAR(x[0], 2.0, 1e-8);
    EXPECT_NEAR(x[1], -3.0, 1e-8);
    EXPECT_NEAR(x[2], 5.0, 1e-8);
    EXPECT_NEAR(finalEnergy, 0.0, 1e-8);
}