#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <algorithm>

/**
* Helper class to build a sparse matrix with a fixed sparsity pattern.
* The matrix is not owned by this class, and must be provided using setMatrix().
* On the first build, the matrix is filled by creating a list of triplets.
* In subsequent builds, entries are inserted into the existing sparse matrix structure.
*/
class FixedPatternMatrixBuilder
{
public:
    FixedPatternMatrixBuilder()
        : m_matrix(nullptr), m_building(false), m_useTriplets(true)
    {
    }

    void setMatrix(Eigen::SparseMatrix<double>* matrix)
    {
        m_matrix = matrix;

        m_building = false;
        m_useTriplets = true;
        m_triplets.clear();
    }

    const Eigen::SparseMatrix<double>* matrix() const
    {
        return m_matrix;
    }

    void begin()
    {
        m_building = true;
        if (!m_useTriplets)
        {
            std::fill_n(m_matrix->valuePtr(), m_matrix->nonZeros(), 0.0);
        }
        else
        {
            m_triplets.clear();
        }
    }

    void end()
    {
        if (m_useTriplets)
        {
            m_matrix->setFromTriplets(m_triplets.begin(), m_triplets.end());
            m_triplets.clear();
            m_triplets.shrink_to_fit();
            m_useTriplets = false;
        }
        m_building = false;
    }

    void addEntry(Eigen::Index row, Eigen::Index col, double v)
    {
#ifndef NDEBUG
        if (!m_building)
        {
            throw std::runtime_error("FixedPatternMatrixBuilder: must call begin() before adding entries");
        }
#endif
        if (m_useTriplets)
        {
            m_triplets.emplace_back(row, col, v);
        }
        else
        {
            m_matrix->coeffRef(row, col) += v;
        }
    }

    template <typename Derived>
    void addBlock(Eigen::Index startRow, Eigen::Index startCol, const Eigen::DenseBase<Derived>& block)
    {
#ifndef NDEBUG
        if (!m_building)
        {
            throw std::runtime_error("FixedPatternMatrixBuilder: must call begin() before adding entries");
        }
#endif
        for (Eigen::Index col = 0; col < block.cols(); ++col)
        {
            const Eigen::Index col0 = startCol + col;
            if (m_useTriplets)
            {
                for (Eigen::Index row = 0; row < block.rows(); ++row)
                {
                    const Eigen::Index row0 = startRow + row;
                    m_triplets.emplace_back(row0, col0, block(row, col));
                }
            }
            else
            {
                const Eigen::Index start = m_matrix->outerIndexPtr()[col0];
                const Eigen::Index end = m_matrix->outerIndexPtr()[col0 + 1];
                const Eigen::Index idx = start == end ? end : m_matrix->data().searchLowerIndex(start, end, startRow);
                for (Eigen::Index row = 0; row < block.rows(); ++row)
                {
                    m_matrix->data().value(idx + row) += block(row, col);
                }
            }
        }
    }

    class ColumnInserter
    {
    public:
        // Row indices must be supplied in ascending order.
        inline void addAscending(Eigen::Index row, double v)
        {
            if (m_builder->m_useTriplets)
            {
                m_builder->m_triplets.emplace_back(row, m_col, v);
                return;
            }

            while (m_innerIndexPtr[m_cursor] < row)
            {
                ++m_cursor;
            }
            m_valuePtr[m_cursor] += v;
        }

    private:
        friend class FixedPatternMatrixBuilder;

        FixedPatternMatrixBuilder* m_builder = nullptr;
        Eigen::Index m_col = 0;
        Eigen::Index m_cursor = 0;
        int* m_innerIndexPtr = nullptr;
        double* m_valuePtr = nullptr;
    };

    ColumnInserter columnInserter(Eigen::Index col)
    {
#ifndef NDEBUG
        if (!m_building)
        {
            throw std::runtime_error("FixedPatternMatrixBuilder: must call begin() before adding entries");
        }
#endif
        ColumnInserter c;
        c.m_builder = this;
        c.m_col = col;
        if (!m_useTriplets)
        {
            c.m_cursor = m_matrix->outerIndexPtr()[col];
            c.m_innerIndexPtr = m_matrix->innerIndexPtr();
            c.m_valuePtr = m_matrix->valuePtr();
        }
        return c;
    }

private:
    Eigen::SparseMatrix<double>* m_matrix;

    bool m_building;
    bool m_useTriplets;
    std::vector<Eigen::Triplet<double>> m_triplets;
};
