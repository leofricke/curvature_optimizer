#include <gtest/gtest.h>

#include "opt/GaussNewtonSystem.h"

TEST(GaussNewtonSystemTest, SystemTest)
{
    GaussNewtonSystem system(2, 2);
    for (size_t i = 0; i < 5; ++i)
    {
        SparseVector row;
        system.beginFilling();

        row.clear();
        system.m_residual[0] = 10.0;
        row.addEntry(0, 1.0);
        system.insertJacobianRow(0, row);

        row.clear();
        system.m_residual[1] = 15.0;
        row.addEntry(0, 2.0);
        row.addEntry(1, 1.0);
        row.addEntry(1, 1.0);
        system.insertJacobianRow(1, row);

        system.endFilling();

        Eigen::Matrix2d expectedMat;
        expectedMat << 5.0, 0.0, 4.0, 4.0;
        Eigen::Vector2d expectedGrad;
        expectedGrad << 40.0, 30.0;

        EXPECT_TRUE(system.m_hessian.isApprox(expectedMat.sparseView()));
        EXPECT_TRUE(system.m_gradient.isApprox(expectedGrad));
    }
}
