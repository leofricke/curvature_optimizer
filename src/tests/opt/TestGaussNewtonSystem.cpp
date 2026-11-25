#include <gtest/gtest.h>

#include "opt/GaussNewtonSystem.h"

TEST(GaussNewtonSystemTest, SystemTest)
{
    GaussNewtonSystem system(2, 2);
    for (size_t i = 0; i < 5; ++i)
    {
        spmat::SpVec<double> row = system.createJacobianRow();
        system.beginFill();

        row.clearPattern();
        row.beginFill();
        system.m_residual[0] = 10.0;
        row.insert(0, 1.0);
        row.endFill();
        system.insertJacobianRow(0, row);

        row.clearPattern();
        row.beginFill();
        system.m_residual[1] = 15.0;
        row.insert(0, 2.0);
        row.insert(1, 1.0);
        row.insert(1, 1.0);
        row.endFill();
        system.insertJacobianRow(1, row);

        const bool patternChanged = system.endFill();
        EXPECT_EQ(patternChanged, i == 0);

        Eigen::Matrix2d expectedMat;
        expectedMat << 5.0, 0.0, 4.0, 4.0;
        Eigen::Vector2d expectedGrad;
        expectedGrad << 40.0, 30.0;

        EXPECT_TRUE(system.m_hessian.map().isApprox(expectedMat.sparseView()));
        EXPECT_TRUE(system.m_gradient.isApprox(expectedGrad));
    }
}

TEST(GaussNewtonSystemTest, PatternGrowsWhenNewEntriesAreInserted)
{
    GaussNewtonSystem system(1, 3);
    spmat::SpVec<double> row = system.createJacobianRow();
    system.m_residual[0] = 1.0;

    system.beginFill();
    row.clearPattern();
    row.beginFill();
    row.insert(0, 1.0);
    row.endFill();
    system.insertJacobianRow(0, row);
    EXPECT_TRUE(system.endFill());

    // same entries: pattern stays the same
    system.beginFill();
    system.insertJacobianRow(0, row);
    EXPECT_FALSE(system.endFill());

    // a new coupling between variables 0 and 2 extends the pattern instead of corrupting the matrix
    system.beginFill();
    row.insert(2, 3.0);
    system.insertJacobianRow(0, row);
    EXPECT_TRUE(system.endFill());

    Eigen::Matrix3d expectedMat;
    expectedMat << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 9.0;
    EXPECT_TRUE(system.m_hessian.map().isApprox(expectedMat.sparseView()));
    EXPECT_TRUE(system.m_hessianPatternDirty);
}
