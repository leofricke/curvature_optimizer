#include <gtest/gtest.h>

#include "opt/GaussNewtonLambdaModule.h"

TEST(GaussNewtonEnergyFunctionTest, SystemTest)
{
    auto module = std::make_shared<GaussNewtonLambdaModule>(
        "Test", 2, 1, [](size_t /*element*/) { return size_t(2); },
        [](size_t /*element*/, GaussNewtonElementContext& ctx)
        {
            ctx.residual(0) = 10.0;
            ctx.residual(1) = 15.0;
            if (ctx.withJacobian())
            {
                ctx.jacobianRow(0).insert(0, 1.0);
                ctx.jacobianRow(1).insert(0, 2.0);
                ctx.jacobianRow(1).insert(1, 1.0);
                ctx.jacobianRow(1).insert(1, 1.0);
            }
        });
    GaussNewtonEnergyFunction energyFunc({ module });
    const Eigen::VectorXd x = Eigen::VectorXd::Zero(2);

    for (size_t i = 0; i < 5; ++i)
    {
        const double energy = energyFunc.evaluateEnergyAndGradient(x);
        EXPECT_EQ(energyFunc.m_hessianPatternDirty, i == 0);

        Eigen::Matrix2d expectedMat;
        expectedMat << 5.0, 0.0, 4.0, 4.0;
        Eigen::Vector2d expectedGrad;
        expectedGrad << 40.0, 30.0;

        EXPECT_TRUE(energyFunc.m_hessian.map().isApprox(expectedMat.sparseView()));
        EXPECT_TRUE(energyFunc.m_gradient.isApprox(expectedGrad));
        EXPECT_DOUBLE_EQ(energy, 0.5 * (10.0 * 10.0 + 15.0 * 15.0));
    }
}

TEST(GaussNewtonEnergyFunctionTest, PatternGrowsWhenNewEntriesAreInserted)
{
    bool coupleVariable2 = false;
    auto module = std::make_shared<GaussNewtonLambdaModule>(
        "Test", 3, 1, [](size_t /*element*/) { return size_t(1); },
        [&coupleVariable2](size_t /*element*/, GaussNewtonElementContext& ctx)
        {
            ctx.residual(0) = 1.0;
            if (ctx.withJacobian())
            {
                ctx.jacobianRow(0).insert(0, 1.0);
                if (coupleVariable2)
                {
                    ctx.jacobianRow(0).insert(2, 3.0);
                }
            }
        });
    GaussNewtonEnergyFunction energyFunc({ module });
    const Eigen::VectorXd x = Eigen::VectorXd::Zero(3);

    energyFunc.evaluateEnergyAndGradient(x);
    EXPECT_TRUE(energyFunc.m_hessianPatternDirty);

    energyFunc.evaluateEnergyAndGradient(x);
    EXPECT_FALSE(energyFunc.m_hessianPatternDirty);

    // a new coupling extends the pattern instead of corrupting the matrix
    coupleVariable2 = true;
    energyFunc.evaluateEnergyAndGradient(x);
    EXPECT_TRUE(energyFunc.m_hessianPatternDirty);

    Eigen::Matrix3d expectedMat;
    expectedMat << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 9.0;
    EXPECT_TRUE(energyFunc.m_hessian.map().isApprox(expectedMat.sparseView()));
}
