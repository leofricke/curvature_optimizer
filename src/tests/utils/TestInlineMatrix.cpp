#include <gtest/gtest.h>

#include "utils/InlineMatrix.h"
#include <Eigen/Core>
#include <utility>

template <typename Scalar, int InlineRows, int InlineCols>
bool usesInlineStorage(const InlineMatrix<Scalar, InlineRows, InlineCols>& m)
{
    const auto* dataPtr = reinterpret_cast<const unsigned char*>(m.data());
    const auto* objStart = reinterpret_cast<const unsigned char*>(&m);
    const auto* objEnd = objStart + sizeof(InlineMatrix<Scalar, InlineRows, InlineCols>);
    return dataPtr >= objStart && dataPtr < objEnd;
}

TEST(InlineMatrixTest, BasicTest)
{
    InlineMatrix<double, 4, 4> m;
    EXPECT_TRUE(m.empty());
    EXPECT_EQ(m.rows(), 0);
    EXPECT_EQ(m.cols(), 0);
    EXPECT_EQ(m.size(), 0);

    m.resize(3, 4);
    EXPECT_FALSE(m.empty());
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 4);
    EXPECT_EQ(m.size(), 12);

    m(0, 0) = 1.0;
    m(1, 2) = 5.5;
    m(2, 3) = 9.9;
    EXPECT_EQ(m(0, 0), 1.0);
    EXPECT_EQ(m(1, 2), 5.5);
    EXPECT_EQ(m(2, 3), 9.9);

    m(1, 2) = 42.0;
    EXPECT_EQ(m(1, 2), 42.0);
}

TEST(InlineMatrixTest, Resize)
{
    InlineMatrix<double, 4, 4> m;

    m.resize(4, 4);
    EXPECT_EQ(m.rows(), 4);
    EXPECT_EQ(m.cols(), 4);
    EXPECT_TRUE(usesInlineStorage(m));

    m.resize(5, 4);
    EXPECT_EQ(m.rows(), 5);
    EXPECT_EQ(m.cols(), 4);
    EXPECT_FALSE(usesInlineStorage(m));

    m.resize(3, 3);
    EXPECT_EQ(m.rows(), 3);
    EXPECT_EQ(m.cols(), 3);
    EXPECT_TRUE(usesInlineStorage(m));

    m.resize(10, 10);
    EXPECT_EQ(m.rows(), 10);
    EXPECT_EQ(m.cols(), 10);
    EXPECT_FALSE(usesInlineStorage(m));

    m.resize(1, 1);
    EXPECT_EQ(m.rows(), 1);
    EXPECT_EQ(m.cols(), 1);
    EXPECT_TRUE(usesInlineStorage(m));
}

TEST(InlineMatrixTest, CopyIsIndependentDeepCopy)
{
    InlineMatrix<double, 2, 2> a;
    a.resize(5, 5);
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            a(i, j) = static_cast<double>(i * 5 + j);
        }
    }

    InlineMatrix<double, 2, 2> b(a);
    InlineMatrix<double, 2, 2> c;
    c = a;

    ASSERT_EQ(b.rows(), 5);
    ASSERT_EQ(b.cols(), 5);
    ASSERT_EQ(c.rows(), 5);
    ASSERT_EQ(c.cols(), 5);
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            EXPECT_EQ(b(i, j), a(i, j));
            EXPECT_EQ(c(i, j), a(i, j));
        }
    }

    b(0, 0) = 999.0;
    c(1, 1) = 888.0;
    EXPECT_EQ(a(0, 0), 0.0);
    EXPECT_EQ(a(1, 1), 6.0);
}

TEST(InlineMatrixTest, MoveTransfersContents)
{
    InlineMatrix<double, 2, 2> a;
    a.resize(5, 5);
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            a(i, j) = static_cast<double>(i * 5 + j);
        }
    }

    const double* originalData = a.data();

    InlineMatrix<double, 2, 2> b(std::move(a));

    EXPECT_EQ(b.data(), originalData);
    ASSERT_EQ(b.rows(), 5);
    for (int i = 0; i < 5; ++i)
    {
        for (int j = 0; j < 5; ++j)
        {
            EXPECT_EQ(b(i, j), static_cast<double>(i * 5 + j));
        }
    }

    EXPECT_EQ(a.rows(), 0);
    EXPECT_EQ(a.cols(), 0);

    InlineMatrix<double, 4, 4> c;
    c.resize(2, 2);
    c(0, 0) = 1.0;
    c(0, 1) = 2.0;
    c(1, 0) = 3.0;
    c(1, 1) = 4.0;
    EXPECT_TRUE(usesInlineStorage(c));

    InlineMatrix<double, 4, 4> d(std::move(c));
    EXPECT_TRUE(usesInlineStorage(d));
    EXPECT_EQ(d(0, 0), 1.0);
    EXPECT_EQ(d(0, 1), 2.0);
    EXPECT_EQ(d(1, 0), 3.0);
    EXPECT_EQ(d(1, 1), 4.0);
    EXPECT_EQ(c.rows(), 0);
    EXPECT_EQ(c.cols(), 0);
}

TEST(InlineMatrixTest, EigenIntegration)
{
    InlineMatrix<double, 4, 4> m;
    m = Eigen::MatrixXd::Identity(3, 3);

    ASSERT_EQ(m.rows(), 3);
    ASSERT_EQ(m.cols(), 3);
    EXPECT_TRUE(usesInlineStorage(m));
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            EXPECT_EQ(m(i, j), i == j ? 1.0 : 0.0);
        }
    }

    const Eigen::Vector3d v { 1.0, 2.0, 3.0 };
    const Eigen::Vector3d result = m.asMap() * v;
    EXPECT_DOUBLE_EQ(result(0), 1.0);
    EXPECT_DOUBLE_EQ(result(1), 2.0);
    EXPECT_DOUBLE_EQ(result(2), 3.0);

    auto sumEntries = [](Eigen::Ref<const Eigen::MatrixXd> mat)
    {
        return mat.sum();
    };
    EXPECT_DOUBLE_EQ(sumEntries(m), 3.0); // trace of identity

    InlineMatrix<double, 4, 4> n;
    n = (m.asMap() * 2.0).eval();
    ASSERT_EQ(n.rows(), 3);
    for (int i = 0; i < 3; ++i)
    {
        EXPECT_DOUBLE_EQ(n(i, i), 2.0);
    }
}

TEST(InlineMatrixTest, ConservativeResize)
{
    InlineMatrix<double, 4, 4> m;

    auto pattern = [](Eigen::Index i, Eigen::Index j)
    {
        return static_cast<double>(i * 100 + j);
    };
    auto fillPattern = [&](InlineMatrix<double, 4, 4>& mat)
    {
        for (Eigen::Index j = 0; j < mat.cols(); ++j)
        {
            for (Eigen::Index i = 0; i < mat.rows(); ++i)
            {
                mat(i, j) = pattern(i, j);
            }
        }
    };
    auto verifyBlock = [&](Eigen::Index rows, Eigen::Index cols)
    {
        for (Eigen::Index j = 0; j < cols; ++j)
        {
            for (Eigen::Index i = 0; i < rows; ++i)
            {
                EXPECT_EQ(m(i, j), pattern(i, j));
            }
        }
    };
    auto step = [&](Eigen::Index rows, Eigen::Index cols, bool expectInline)
    {
        const Eigen::Index prevRows = m.rows();
        const Eigen::Index prevCols = m.cols();
        const Eigen::Index keepRows = prevRows < rows ? prevRows : rows;
        const Eigen::Index keepCols = prevCols < cols ? prevCols : cols;

        m.conservativeResize(rows, cols);

        EXPECT_EQ(m.rows(), rows);
        EXPECT_EQ(m.cols(), cols);
        EXPECT_EQ(usesInlineStorage(m), expectInline);
        verifyBlock(keepRows, keepCols);

        for (Eigen::Index j = 0; j < cols; ++j)
        {
            for (Eigen::Index i = 0; i < rows; ++i)
            {
                if (i >= keepRows || j >= keepCols)
                {
                    m(i, j) = -1.0;
                }
            }
        }
        verifyBlock(keepRows, keepCols);

        fillPattern(m);
    };

    m.resize(3, 3);
    fillPattern(m);

    const double* dataBefore = m.data();
    m.conservativeResize(3, 3);
    EXPECT_EQ(m.data(), dataBefore);
    EXPECT_TRUE(usesInlineStorage(m));
    verifyBlock(3, 3);

    step(3, 4, true);
    step(3, 2, true);
    step(4, 2, true);
    step(2, 4, true);
    step(5, 5, false);
    step(6, 3, false);
    step(3, 3, true);
    step(1, 1, true);
}