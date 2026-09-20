#pragma once

#include "utils/InlineVector.h"
#include "utils/InlineMatrix.h"
#include <pmp/surface_mesh.h>

template <typename GradientType, int ExpectedSize = 4>
    requires(std::is_base_of<Eigen::MatrixBase<GradientType>, GradientType>::value)
class LocalGradient
{
    static_assert(GradientType::RowsAtCompileTime > 0, "GradientType must have compile-time known number of rows");
    static_assert(GradientType::ColsAtCompileTime > 0, "GradientType must have compile-time known number of cols");
    static_assert(ExpectedSize > 1, "ExpectedSize must be greater than 1");

public:
    using VertexVectorType = InlineVector<pmp::Vertex, ExpectedSize>;
    using GradientMatrixType = InlineMatrix<double, GradientType::RowsAtCompileTime, ExpectedSize * GradientType::ColsAtCompileTime>;
    using BlockType = Eigen::Block<typename GradientMatrixType::MapType, GradientType::RowsAtCompileTime, GradientType::ColsAtCompileTime>;
    using ConstBlockType = Eigen::Block<typename GradientMatrixType::ConstMapType, GradientType::RowsAtCompileTime, GradientType::ColsAtCompileTime>;

    LocalGradient()
    {
        m_gradient.resize(GradientType::RowsAtCompileTime, ExpectedSize * GradientType::ColsAtCompileTime);
        clear();
    }

    void clear()
    {
        m_vertices.clear();
        m_gradient.asMap().setZero();
    }

    int addVertex(pmp::Vertex v)
    {
        m_vertices.emplace_back(v);
        if (m_gradient.cols() < static_cast<Eigen::Index>(m_vertices.size()) * GradientType::ColsAtCompileTime)
        {
            const Eigen::Index prevRows = m_gradient.rows();
            const Eigen::Index prevCols = m_gradient.cols();
            m_gradient.conservativeResize(prevRows, 2 * prevCols);
            m_gradient.asMap().block(0, prevCols, prevRows, prevCols).setZero();
        }
        return m_vertices.size() - 1;
    }

    int findVertexIdx(pmp::Vertex v) const
    {
        for (size_t idx = 0; idx < m_vertices.size(); ++idx)
        {
            if (m_vertices[idx] == v)
            {
                return idx;
            }
        }
        return -1;
    }

    bool containsVertex(pmp::Vertex v) const
    {
        return findVertexIdx(v) >= 0;
    }

    BlockType blockByIdx(int idx)
    {
        return m_gradient.asMap().template block<GradientType::RowsAtCompileTime, GradientType::ColsAtCompileTime>(0, idx * GradientType::ColsAtCompileTime);
    }

    ConstBlockType blockByIdx(int idx) const
    {
        return m_gradient.asMap().template block<GradientType::RowsAtCompileTime, GradientType::ColsAtCompileTime>(0, idx * GradientType::ColsAtCompileTime);
    }

    BlockType operator[](pmp::Vertex v)
    {
        return this->operator()(v);
    }

    ConstBlockType operator[](pmp::Vertex v) const
    {
        return this->operator()(v);
    }

    BlockType operator()(pmp::Vertex v)
    {
        int idx = findVertexIdx(v);
        if (idx < 0)
        {
            idx = addVertex(v);
        }
        return blockByIdx(idx);
    }

    ConstBlockType operator()(pmp::Vertex v) const
    {
        const int idx = findVertexIdx(v);
        if (idx < 0)
        {
            throw std::runtime_error("Vertex does not exists");
        }
        return blockByIdx(idx);
    }

    VertexVectorType::const_iterator begin() const
    {
        return m_vertices.begin();
    }

    VertexVectorType::const_iterator end() const
    {
        return m_vertices.end();
    }

    LocalGradient& operator*=(double s)
    {
        m_gradient.asMap() *= s;
        return *this;
    }

    LocalGradient& operator/=(double s)
    {
        m_gradient.asMap() /= s;
        return *this;
    }

    template <int N>
    LocalGradient& multiplyAdd(double a, const LocalGradient<GradientType, N>& rhs)
    {
        for (pmp::Vertex v : rhs)
        {
            int idx = findVertexIdx(v);
            if (idx < 0)
            {
                idx = addVertex(v);
            }
            blockByIdx(idx) += a * rhs[v];
        }
        return *this;
    }

    template <int N>
    LocalGradient& operator+=(const LocalGradient<GradientType, N>& rhs)
    {
        return multiplyAdd(1.0, rhs);
    }

    template <int N>
    LocalGradient& operator-=(const LocalGradient<GradientType, N>& rhs)
    {
        return multiplyAdd(-1.0, rhs);
    }

private:
    VertexVectorType m_vertices;
    GradientMatrixType m_gradient;
};
