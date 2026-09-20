#pragma once

#include "LocalGradient.h"

template <typename T>
concept EigenMatrixType = std::is_base_of<Eigen::MatrixBase<T>, T>::value;

template <EigenMatrixType T, int N, int M, typename DerivFunc>
void applyChainRule1(double x, DerivFunc&& d1, const LocalGradient<T, N>& gx, LocalGradient<T, M>& gf)
{
    const double d1x = d1(x);
    gf.multiplyAdd(d1x, gx);
}

template <EigenMatrixType T, int N, typename DerivFunc>
void applyChainRule1Inplace(double x, DerivFunc&& d1, LocalGradient<T, N>& gx)
{
    const double d1x = d1(x);
    gx *= d1x;
}

template <EigenMatrixType T, int N, int M, int O>
void applyProductRule1(double x, double y, const LocalGradient<T, N>& gx, const LocalGradient<T, M>& gy, LocalGradient<T, O>& gf)
{
    // first order gradient terms
    gf.multiplyAdd(y, gx);
    gf.multiplyAdd(x, gy);
}

template <EigenMatrixType T, int N, int M>
void applyProductRule1Inplace(double x, double y, const LocalGradient<T, N>& gy, LocalGradient<T, M>& gx)
{
    // first order gradient terms
    gx *= y;
    gx.multiplyAdd(x, gy);
}

template <EigenMatrixType T, int N, int M, int O>
void applyQuotientRule1(double x, double y, const LocalGradient<T, N>& gx, const LocalGradient<T, M>& gy, LocalGradient<T, O>& gf)
{
    const double inv = 1.0 / y;
    const double inv2 = 1.0 / (y * y);
    gf.multiplyAdd(inv, gx);
    gf.multiplyAdd(-x * inv2, gy);
}

template <EigenMatrixType T, int N, int M>
void applyQuotientRule1Inplace(double x, double y, const LocalGradient<T, N>& gy, LocalGradient<T, M>& gx)
{
    const double inv = 1.0 / y;
    gx.multiplyAdd(-x * inv, gy);
    gx *= inv;
}
