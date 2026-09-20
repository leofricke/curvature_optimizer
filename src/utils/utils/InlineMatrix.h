#pragma once

#include <Eigen/Core>
#include <new>
#include <cassert>
#include <algorithm>

template <typename Scalar, int InlineRows, int InlineCols>
class InlineMatrix
{
    static_assert(InlineRows >= 1 && InlineCols >= 1, "InlineMatrix requires inline dimensions of at least 1.");

public:
    static constexpr std::size_t InlineCapacity = static_cast<std::size_t>(InlineRows) * static_cast<std::size_t>(InlineCols);

    using MapType = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>;
    using ConstMapType = Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>;

    InlineMatrix() noexcept
        : m_data(m_inline), m_rows(0), m_cols(0)
    {
    }

    ~InlineMatrix() noexcept
    {
        if (!using_inline())
        {
            free_heap();
        }
    }

    InlineMatrix(const InlineMatrix& other) noexcept
        : m_data(m_inline), m_rows(0), m_cols(0)
    {
        copy_from(other);
    }

    InlineMatrix& operator=(const InlineMatrix& other) noexcept
    {
        if (this != &other)
        {
            copy_from(other);
        }
        return *this;
    }

    InlineMatrix(InlineMatrix&& other) noexcept
        : m_data(m_inline), m_rows(0), m_cols(0)
    {
        move_from(std::move(other));
    }

    InlineMatrix& operator=(InlineMatrix&& other) noexcept
    {
        if (this != &other)
        {
            if (!using_inline())
            {
                free_heap();
            }
            move_from(std::move(other));
        }
        return *this;
    }

    template <typename OtherDerived>
    explicit InlineMatrix(const Eigen::MatrixBase<OtherDerived>& other)
        : m_data(m_inline), m_rows(0), m_cols(0)
    {
        resize(other.rows(), other.cols());
        asMap() = other;
    }

    template <typename OtherDerived>
    InlineMatrix& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        resize(other.rows(), other.cols());
        asMap() = other;
        return *this;
    }

    void resize(Eigen::Index rows, Eigen::Index cols)
    {
        assert(rows >= 0 && cols >= 0);
        const std::size_t needed = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);

        if (needed <= InlineCapacity)
        {
            if (!using_inline())
            {
                free_heap();
                m_data = m_inline;
            }
        }
        else
        {
            if (!using_inline())
            {
                free_heap();
            }
            allocate_heap(needed);
        }

        m_rows = rows;
        m_cols = cols;
    }

    void conservativeResize(Eigen::Index rows, Eigen::Index cols)
    {
        assert(rows >= 0 && cols >= 0);
        if (rows == m_rows && cols == m_cols)
        {
            return;
        }
        const std::size_t needed = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
        const bool needsHeap = needed > InlineCapacity;
        if (rows == m_rows && needed <= InlineCapacity && using_inline())
        {
            m_cols = cols;
            return;
        }
        const Eigen::Index minRows = std::min(rows, m_rows);
        const Eigen::Index minCols = std::min(cols, m_cols);

        if (!needsHeap && using_inline())
        {
            Scalar temp[InlineCapacity];
            copy_block(temp, rows, m_data, m_rows, minRows, minCols);
            copy_block(m_inline, rows, temp, rows, minRows, minCols);
        }
        else
        {
            Scalar* dst = needsHeap ? static_cast<Scalar*>(::operator new(sizeof(Scalar) * needed, std::align_val_t { EIGEN_DEFAULT_ALIGN_BYTES })) : m_inline;
            copy_block(dst, rows, m_data, m_rows, minRows, minCols);
            if (!using_inline())
            {
                free_heap();
            }
            m_data = dst;
        }

        m_rows = rows;
        m_cols = cols;
    }

    // clang-format off
    Eigen::Index rows()  const { return m_rows; }
    Eigen::Index cols()  const { return m_cols; }
    Eigen::Index size()  const { return m_rows * m_cols; }
    bool empty() const { return m_rows == 0 || m_cols == 0; }
    // clang-format on

    Scalar& operator()(Eigen::Index i, Eigen::Index j)
    {
        assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
        return m_data[flat_index(i, j)];
    }

    const Scalar& operator()(Eigen::Index i, Eigen::Index j) const
    {
        assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
        return m_data[flat_index(i, j)];
    }

    Scalar* data()
    {
        return m_data;
    }

    const Scalar* data() const
    {
        return m_data;
    }

    MapType asMap()
    {
        return MapType(m_data, m_rows, m_cols);
    }

    ConstMapType asMap() const
    {
        return ConstMapType(m_data, m_rows, m_cols);
    }

    operator Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>()
    {
        return Eigen::Ref<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(asMap());
    }

    operator Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>() const
    {
        return Eigen::Ref<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>(asMap());
    }

private:
    alignas(EIGEN_DEFAULT_ALIGN_BYTES) Scalar m_inline[InlineCapacity];
    Scalar* m_data;
    Eigen::Index m_rows;
    Eigen::Index m_cols;

    bool using_inline() const
    {
        return m_data == m_inline;
    }

    void allocate_heap(std::size_t count)
    {
        m_data = static_cast<Scalar*>(::operator new(sizeof(Scalar) * count, std::align_val_t { EIGEN_DEFAULT_ALIGN_BYTES }));
    }

    void free_heap()
    {
        ::operator delete(m_data, std::align_val_t { EIGEN_DEFAULT_ALIGN_BYTES });
    }

    Eigen::Index flat_index(Eigen::Index i, Eigen::Index j) const
    {
        return j * m_rows + i; // column-major
    }

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
    void copy_from(const InlineMatrix& other)
    {
        m_rows = other.m_rows;
        m_cols = other.m_cols;

        if (!using_inline())
        {
            free_heap();
        }
        if (other.using_inline())
        {
            m_data = m_inline;
        }
        else
        {
            allocate_heap(static_cast<std::size_t>(m_rows) * static_cast<std::size_t>(m_cols));
        }

        std::copy_n(other.m_data, static_cast<std::size_t>(m_rows) * static_cast<std::size_t>(m_cols), m_data);
    }
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

    void move_from(InlineMatrix&& other)
    {
        m_rows = other.m_rows;
        m_cols = other.m_cols;

        if (other.using_inline())
        {
            m_data = m_inline;
            std::copy_n(other.m_data, static_cast<std::size_t>(m_rows) * static_cast<std::size_t>(m_cols), m_data);
        }
        else
        {
            m_data = other.m_data;
            other.m_data = other.m_inline;
        }

        other.m_rows = 0;
        other.m_cols = 0;
    }

    static void copy_block(Scalar* dst, Eigen::Index dstStride, const Scalar* src, Eigen::Index srcStride, Eigen::Index minRows, Eigen::Index minCols)
    {
        if (srcStride == dstStride && minRows == dstStride)
        {
            std::copy_n(src, static_cast<std::size_t>(dstStride) * static_cast<std::size_t>(minCols), dst);
        }
        else
        {
            for (Eigen::Index j = 0; j < minCols; ++j)
            {
                std::copy_n(src + j * srcStride, minRows, dst + j * dstStride);
            }
        }
    }
};