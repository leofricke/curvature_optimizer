#pragma once

#include <algorithm>
#include <cassert>
#include <limits>
#include <spmat/SpMat.h>

template <typename Scalar, spmat::SpMatLayout Layout, spmat::SpMatSymmetry Symmetry, spmat::SpMatThreadSafety ThreadSafety, typename Index>
Scalar minDiagonal(const spmat::SpMat<Scalar, Layout, Symmetry, ThreadSafety, Index>& mat)
{
    Scalar min = std::numeric_limits<Scalar>::max();
    for (Index i = 0; i < mat.rows(); ++i)
    {
        min = std::min(min, mat.diagonal(i));
    }
    return min;
}

template <typename Scalar, spmat::SpMatLayout Layout, spmat::SpMatSymmetry Symmetry, spmat::SpMatThreadSafety ThreadSafety, typename Index>
void regularize(spmat::SpMat<Scalar, Layout, Symmetry, ThreadSafety, Index>& mat, Scalar regularizer)
{
    if (regularizer > Scalar(0))
    {
        for (Index i = 0; i < mat.rows(); ++i)
        {
            Scalar* d = mat.diagonalPtr(i);
            assert(d && "regularize() requires the diagonal to be part of the sparsity pattern");
            *d += regularizer;
        }
    }
}
