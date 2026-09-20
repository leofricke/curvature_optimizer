#include <gtest/gtest.h>

#include "utils/KoppEigendecomposition.h"

#include <Eigen/Eigenvalues>
#include <algorithm>

inline Eigen::Vector3d sorted(const Eigen::Vector3d& v)
{
    Eigen::Vector3d s = v;
    std::sort(s.data(), s.data() + 3);
    return s;
}

inline bool isSorted(const Eigen::Vector3d& v)
{
    return std::is_sorted(v.data(), v.data() + 3);
}

inline bool approxEqualUnordered(const Eigen::Vector3d& a, const Eigen::Vector3d& b, double tol = 1e-10)
{
    return (sorted(a) - sorted(b)).norm() < tol;
}

TEST(KoppEigendecomposition, DiagonalMatrixEigenvalues)
{
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
    A(0, 0) = 10.0;
    A(1, 1) = 2.0;
    A(2, 2) = 30.0;

    Eigen::Vector3d w;
    eigendecomposition3(A, w);

    EXPECT_TRUE(isSorted(w));

    Eigen::Vector3d w_ref(10.0, 2.0, 30.0);
    EXPECT_TRUE(approxEqualUnordered(w, w_ref));
}

TEST(KoppEigendecomposition, DiagonalMatrixEigenvectors)
{
    Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
    A(0, 0) = 4.0;
    A(1, 1) = 5.0;
    A(2, 2) = 6.0;

    Eigen::Vector3d w;
    Eigen::Matrix3d Q;
    eigendecomposition3(A, w, Q);

    EXPECT_TRUE(isSorted(w));

    Eigen::Vector3d w_ref(4.0, 5.0, 6.0);
    EXPECT_TRUE(approxEqualUnordered(w, w_ref));

    Eigen::Matrix3d A_rec = Q * w.asDiagonal() * Q.transpose();
    EXPECT_TRUE(A_rec.isApprox(A, 1e-12));
}

TEST(KoppEigendecomposition, RandomSymmetricMatrix)
{
    Eigen::Matrix3d M = Eigen::Matrix3d::Random();
    Eigen::Matrix3d A = 0.5 * (M + M.transpose());

    Eigen::Vector3d w;
    Eigen::Matrix3d Q;
    eigendecomposition3(A, w, Q);

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(A);
    ASSERT_EQ(es.info(), Eigen::Success);

    EXPECT_TRUE(approxEqualUnordered(w, es.eigenvalues()));

    Eigen::Matrix3d A_rec = Q * w.asDiagonal() * Q.transpose();
    EXPECT_TRUE(A_rec.isApprox(A, 1e-12));
}

TEST(KoppEigendecomposition, RankDeficientMatrix)
{
    Eigen::Matrix3d A;
    A << 1, 0, 0, 0, 2, 0, 0, 0, 0;

    Eigen::Vector3d w;
    Eigen::Matrix3d Q;
    eigendecomposition3(A, w, Q);

    Eigen::Vector3d w_ref(0.0, 1.0, 2.0);
    EXPECT_TRUE(approxEqualUnordered(w, w_ref));

    Eigen::Matrix3d A_rec = Q * w.asDiagonal() * Q.transpose();
    EXPECT_TRUE(A_rec.isApprox(A, 1e-12));
}

TEST(KoppEigendecomposition, EigenvectorsSatisfyAvEqualsLambdaV)
{
    Eigen::Matrix3d A;
    A << 2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0;

    Eigen::Vector3d w;
    Eigen::Matrix3d Q;
    eigendecomposition3(A, w, Q);

    for (int i = 0; i < 3; ++i)
    {
        Eigen::Vector3d v = Q.col(i);
        EXPECT_TRUE((A * v).isApprox(w(i) * v, 1e-12));
    }
}
