#pragma once

#include <Eigen/Sparse>

template <typename SparseMat>
typename SparseMat::Scalar minDiagonal(const SparseMat& mat)
{
    auto min = std::numeric_limits<typename SparseMat::Scalar>::max();
    for (Eigen::Index i = 0; i < mat.cols(); ++i)
    {
        const auto d = mat.coeff(i, i);
        if (d < min)
        {
            min = d;
        }
    }
    return min;
}

template <typename SparseMat>
void regularize(SparseMat&& mat, double regularizer)
{
    if (regularizer > 0.0)
    {
        for (Eigen::Index i = 0; i < mat.cols(); ++i)
        {
            mat.coeffRef(i, i) += regularizer;
        }
    }
}

template <typename SparseMat>
void regularize(SparseMat&& mat, const Eigen::Vector<typename SparseMat::Scalar, Eigen::Dynamic>& diag)
{
    if (mat.rows() == diag.size() && mat.cols() == diag.size())
    {
        for (Eigen::Index i = 0; i < mat.cols(); ++i)
        {
            mat.coeffRef(i, i) += diag[i];
        }
    }
}

template <typename Scalar>
Eigen::Matrix<Scalar, 3, 3> crossMatrix(const Eigen::Vector<Scalar, 3>& a)
{
    Eigen::Matrix<Scalar, 3, 3> m;
    // clang-format off
    m << 0.0, -a[2], a[1],
         a[2], 0.0, -a[0],
        -a[1], a[0], 0.0;
    // clang-format on
    return m;
}

template <typename Derived>
typename Derived::PlainObject sym(const Eigen::MatrixBase<Derived>& m)
{
    typename Derived::PlainObject p = m;
    return (p + p.transpose()) / 2.0;
}

template <typename Scalar>
void computeOrthoBasis(const Eigen::Vector<Scalar, 3>& v, Eigen::Vector<Scalar, 3>& e0, Eigen::Vector<Scalar, 3>& e1)
{
    // pick a non-parallel vector
    if (std::abs(v[0]) > std::abs(v[1]))
    {
        e0 << -v[2], 0.0, v[0];
    }
    else
    {
        e0 << 0.0, v[2], -v[1];
    }

    // gram-schmidt process
    e0 = (e0 - e0.dot(v) * v).normalized();
    e1 = v.cross(e0).normalized();
}
