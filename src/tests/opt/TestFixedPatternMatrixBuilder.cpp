#include <gtest/gtest.h>

#include "opt/utils/FixedPatternMatrixBuilder.h"

TEST(FixedPatternMatrixBuilderTest, BuilderTest)
{
    Eigen::SparseMatrix<double> matrix(5, 5);
    FixedPatternMatrixBuilder builder;
    builder.setMatrix(&matrix);

    for (size_t i = 0; i < 5; ++i)
    {
        builder.begin();
        builder.addEntry(1, 0, 1.0);
        builder.addEntry(2, 3, 2.0);
        builder.addEntry(4, 2, 3.0);
        const Eigen::Matrix2d m = Eigen::Matrix2d::Ones();
        builder.addBlock(3, 3, m);
        builder.end();

        EXPECT_EQ(matrix.sum(), 10.0);
        EXPECT_EQ(matrix.coeff(1, 0), 1.0);
        EXPECT_EQ(matrix.coeff(2, 3), 2.0);
        EXPECT_EQ(matrix.coeff(4, 2), 3.0);
        EXPECT_EQ(matrix.coeff(3, 3), 1.0);
        EXPECT_EQ(matrix.coeff(4, 3), 1.0);
        EXPECT_EQ(matrix.coeff(3, 4), 1.0);
        EXPECT_EQ(matrix.coeff(4, 4), 1.0);
    }
}
