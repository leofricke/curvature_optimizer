#pragma once

#include <Eigen/Core>

template <typename T, typename R, int N>
using NumericalDiffFunc = std::function<R(const std::array<T, N>&)>;

template <typename T>
concept IsEigenVector = Eigen::internal::traits<std::remove_cv_t<T>>::RowsAtCompileTime != Eigen::Dynamic && Eigen::internal::traits<std::remove_cv_t<T>>::ColsAtCompileTime == 1;

template <typename T, typename R, int N>
    requires IsEigenVector<T>
void numericalDeriv(NumericalDiffFunc<T, R, N> func, const std::array<T, N>& args, std::array<std::array<R, T::RowsAtCompileTime>, N>& derivs, double eps = 1e-8)
{
    std::array<T, N> args_plus;
    std::array<T, N> args_minus;
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < T::RowsAtCompileTime; ++j)
        {
            args_plus = args;
            args_minus = args;
            args_plus[i][j] += eps;
            args_minus[i][j] -= eps;

            derivs[i][j] = (func(args_plus) - func(args_minus)) / (2.0 * eps);
        }
    }
}

template <typename Scalar>
bool checkDeriv(Scalar analytic, Scalar numerical, Scalar atol = 1e-8, Scalar rtol = 1e-7)
{
    const Scalar e = std::abs(analytic - numerical);
    if (e >= atol + rtol * std::abs(numerical))
    {
        std::cout << "analytic: " << analytic << "  numerical: " << numerical << std::endl;
        return false;
    }
    return true;
}

template <typename Derived1, typename Derived2>
bool checkDerivMatrix(const Eigen::DenseBase<Derived1>& analytic, const Eigen::DenseBase<Derived2>& numerical, typename Derived1::Scalar atol = 1e-8, typename Derived1::Scalar rtol = 1e-7)
{
    for (Eigen::Index r = 0; r < analytic.rows(); ++r)
    {
        for (Eigen::Index c = 0; c < analytic.cols(); ++c)
        {
            const typename Derived1::Scalar e = std::abs(analytic(r, c) - numerical(r, c));
            if (e >= atol + rtol * std::abs(numerical(r, c)))
            {
                std::cout << "analytic: " << analytic(r, c) << "  numerical: " << numerical(r, c) << std::endl;
                return false;
            }
        }
    }
    return true;
}