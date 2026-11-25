#pragma once

#include "DefaultInitAllocator.h"
#include "SpVec.h"

#include <Eigen/SparseCore>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <limits>
#include <mutex>
#include <type_traits>
#include <vector>

namespace spmat
{

enum class SpMatLayout
{
    ColMajor,
    RowMajor
};

// Lower/Upper: the matrix is symmetric and only that triangle is stored. All insert operations drop entries outside it,
// so callers may insert any contributions. Coefficient access only covers the stored triangle.
enum class SpMatSymmetry
{
    None,
    Lower,
    Upper
};

// Enabled: insert(), insertBlock() and insertOuterProduct() may be called concurrently from several threads during a fill.
// Disabled: no locks or atomics; a single insertOuterProduct(..., parallel = true) call is still safe.
enum class SpMatThreadSafety
{
    Disabled,
    Enabled
};

template <typename Scalar, SpMatLayout Layout = SpMatLayout::ColMajor, SpMatSymmetry Symmetry = SpMatSymmetry::None, SpMatThreadSafety ThreadSafety = SpMatThreadSafety::Disabled, typename Index = int>
class SpMat
{
    static_assert(std::is_integral_v<Index> && std::is_signed_v<Index>, "IndexType must be a signed integral type");
    static_assert(Layout == SpMatLayout::ColMajor || Layout == SpMatLayout::RowMajor, "Layout must be either ColMajor or RowMajor");

public:
    static constexpr bool IsColMajor = Layout == SpMatLayout::ColMajor;
    static constexpr bool IsRowMajor = Layout == SpMatLayout::RowMajor;
    static constexpr bool IsSymmetric = Symmetry != SpMatSymmetry::None;
    static constexpr bool IsThreadSafe = ThreadSafety == SpMatThreadSafety::Enabled;

    using EigenType = Eigen::SparseMatrix<Scalar, IsColMajor ? Eigen::ColMajor : Eigen::RowMajor, Index>;
    using MapType = Eigen::Map<EigenType>;
    using ConstMapType = Eigen::Map<const EigenType>;

    struct BucketEntry
    {
        Index m_inner;
        Scalar m_value;
    };
    using Bucket = std::vector<BucketEntry>;
    using IndexArray = DefaultInitVector<Index>;
    using ScalarArray = DefaultInitVector<Scalar>;

    SpMat()
        : SpMat(0, 0)
    {
    }

    SpMat(Index rows, Index cols)
    {
        assert(rows >= 0 && cols >= 0);
        assert((!IsSymmetric || rows == cols) && "symmetric SpMat must be square");
        m_rows = rows;
        m_cols = cols;
        m_outer.assign(static_cast<size_t>(outerSize()) + 1, Index(0));
        resetBucketStorage();
    }

    inline void resize(Index rows, Index cols)
    {
        assert(!m_filling && "resize() must not be called between beginFill()/endFill()");
        assert(rows >= 0 && cols >= 0);
        assert((!IsSymmetric || rows == cols) && "symmetric SpMat must be square");
        if (rows == m_rows && cols == m_cols)
        {
            return;
        }

        m_rows = rows;
        m_cols = cols;
        m_inner.clear();
        m_values.clear();
        m_outer.assign(static_cast<size_t>(outerSize()) + 1, Index(0));
        resetBucketStorage();
    }

    inline void clearPattern()
    {
        assert(!m_filling && "clearPattern() must not be called between beginFill()/endFill()");
        m_inner.clear();
        m_values.clear();
        std::fill(m_outer.begin(), m_outer.end(), Index(0));
        for (Bucket& b : m_buckets)
        {
            b.clear();
        }
    }

    inline void setZero()
    {
        std::fill(m_values.begin(), m_values.end(), Scalar(0.0));
    }

    inline Index rows() const
    {
        return m_rows;
    }

    inline Index cols() const
    {
        return m_cols;
    }

    inline Index outerSize() const
    {
        return IsRowMajor ? m_rows : m_cols;
    }

    inline Index innerSize() const
    {
        return IsRowMajor ? m_cols : m_rows;
    }

    inline Index nonZeros() const
    {
        return static_cast<Index>(m_inner.size());
    }

    inline const IndexArray& outerArray() const
    {
        return m_outer;
    }

    inline const IndexArray& innerArray() const
    {
        return m_inner;
    }

    inline const ScalarArray& valuesArray() const
    {
        return m_values;
    }

    inline bool isFilling() const
    {
        return m_filling;
    }

    inline size_t pendingBucketEntries() const
    {
        size_t pending = 0;
        for (const Bucket& bucket : m_buckets)
        {
            pending += bucket.size();
        }
        return pending;
    }

    inline void beginFill()
    {
        assert(!m_filling && "beginFill() called twice without endFill()");

        setZero();
        m_filling = true;
    }

    inline bool endFill()
    {
        assert(m_filling && "endFill() without a matching beginFill()");
        m_filling = false;

        const bool hasPending = std::any_of(m_buckets.begin(), m_buckets.end(), [](const Bucket& b) { return !b.empty(); });
        if (!hasPending)
        {
            return false;
        }
        mergeBuckets();
        return true;
    }

    inline void insert(Index row, Index col, Scalar value)
    {
        assert(m_filling && "insert() must be called between beginFill()/endFill()");
        assert(row >= 0 && row < m_rows && col >= 0 && col < m_cols);

        const Index outer = IsRowMajor ? row : col;
        const Index inner = IsRowMajor ? col : row;
        if ((OuterStartsAtDiagonal && inner < outer) || (OuterEndsAtDiagonal && inner > outer))
        {
            return;
        }
        const Index idx = findIndexOuterInner(outer, inner);

        std::unique_lock<std::mutex> lock = makeOuterLock(outer);
        if (idx >= 0)
        {
            addValue(m_values[static_cast<size_t>(idx)], value);
        }
        else
        {
            ensureOuterLocked(lock, outer);
            m_buckets[static_cast<size_t>(outer)].push_back(BucketEntry { inner, value });
        }
    }

    template <class Derived>
    inline void insertBlock(Index startRow, Index startCol, const Eigen::MatrixBase<Derived>& block)
    {
        assert(m_filling && "insertBlock() must be called between beginFill()/endFill()");
        assert(startRow >= 0 && startRow + block.rows() <= m_rows && startCol >= 0 && startCol + block.cols() <= m_cols);

        const Index blockOuter = static_cast<Index>(IsRowMajor ? block.rows() : block.cols());
        const Index blockInner = static_cast<Index>(IsRowMajor ? block.cols() : block.rows());
        const Index startOuter = IsRowMajor ? startRow : startCol;
        const Index startInner = IsRowMajor ? startCol : startRow;
        const Index* const innerPtr = m_inner.data();
        Scalar* const valuePtr = m_values.data();
        for (Index o = 0; o < blockOuter; ++o)
        {
            const Index outer = startOuter + o;
            const Index innerBegin = OuterStartsAtDiagonal ? std::max(startInner, outer) : startInner;
            const Index innerEnd = OuterEndsAtDiagonal ? std::min(startInner + blockInner, outer + 1) : startInner + blockInner;
            if (innerBegin >= innerEnd)
            {
                continue;
            }
            const Index* const last = innerPtr + m_outer[static_cast<size_t>(outer) + 1];
            const Index* it = std::lower_bound(innerPtr + m_outer[static_cast<size_t>(outer)], last, innerBegin);
            std::unique_lock<std::mutex> lock = makeOuterLock(outer);
            for (Index inner = innerBegin; inner < innerEnd; ++inner)
            {
                const Index i = inner - startInner;
                const Scalar value = IsRowMajor ? block(o, i) : block(i, o);
                while (it != last && *it < inner)
                {
                    ++it;
                }
                if (it != last && *it == inner)
                {
                    addValue(valuePtr[it - innerPtr], value);
                }
                else
                {
                    ensureOuterLocked(lock, outer);
                    m_buckets[static_cast<size_t>(outer)].push_back(BucketEntry { inner, value });
                }
            }
        }
    }

    inline void insertOuterProduct(const SpVec<Scalar, Index>& v, Scalar w = 1.0, bool parallel = false)
    {
        assert(m_filling && "insertOuterProduct() must be called between beginFill()/endFill()");
        assert(!v.isFilling() && "insertOuterProduct() requires a sorted SpVec; call v.endFill() first");

        const auto* const entries = v.entries().data();
        const Index nnz = v.nonZeros();

        const auto insertOuter = [=, this](Index o)
        {
            const Index outer = entries[o].m_idx;
            assert(outer >= 0 && outer < outerSize());
            const Index end = m_outer[static_cast<size_t>(outer) + 1];
            const Index stop = OuterEndsAtDiagonal ? o + 1 : nnz;
            const Scalar rowScale = entries[o].m_value * w;
            const Index* const innerPtr = m_inner.data();
            Scalar* const valuePtr = m_values.data();

            std::unique_lock<std::mutex> lock = makeOuterLock(outer);
            const Index first = OuterStartsAtDiagonal ? o : 0;
            Index p = static_cast<Index>(std::lower_bound(innerPtr + m_outer[static_cast<size_t>(outer)], innerPtr + end, entries[first].m_idx) - innerPtr);
            for (Index i = first; i < stop; ++i)
            {
                const Index inner = entries[i].m_idx;
                assert(inner >= 0 && inner < innerSize());
                const Scalar value = rowScale * entries[i].m_value;

                while (p < end && innerPtr[p] < inner)
                {
                    ++p;
                }
                if (p < end && innerPtr[p] == inner)
                {
                    addValue(valuePtr[p], value);
                }
                else
                {
                    ensureOuterLocked(lock, outer);
                    m_buckets[static_cast<size_t>(outer)].push_back(BucketEntry { inner, value });
                }
            }
        };

        if (parallel)
        {
#pragma omp parallel for schedule(guided)
            for (Index o = 0; o < nnz; ++o)
            {
                insertOuter(o);
            }
        }
        else
        {
            for (Index o = 0; o < nnz; ++o)
            {
                insertOuter(o);
            }
        }
    }

    inline Scalar coeff(Index row, Index col) const
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        const Index idx = findIndex(row, col);
        return idx >= 0 ? m_values[idx] : Scalar(0.0);
    }

    inline Scalar& coeffRef(Index row, Index col)
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        const Index idx = findIndex(row, col);
        assert(idx >= 0 && "coeffRef() on a coefficient that is not in the pattern");
        return m_values[idx];
    }

    inline Scalar* coeffPtr(Index row, Index col)
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        const Index idx = findIndex(row, col);
        return idx >= 0 ? m_values.data() + idx : nullptr;
    }

    inline const Scalar* coeffPtr(Index row, Index col) const
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        const Index idx = findIndex(row, col);
        return idx >= 0 ? m_values.data() + idx : nullptr;
    }

    inline Scalar diagonal(Index i) const
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        const Index idx = findDiagonalIndex(i);
        return idx >= 0 ? m_values[static_cast<size_t>(idx)] : Scalar(0.0);
    }

    inline Scalar* diagonalPtr(Index i)
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        const Index idx = findDiagonalIndex(i);
        return idx >= 0 ? m_values.data() + idx : nullptr;
    }

    inline const Scalar* diagonalPtr(Index i) const
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        const Index idx = findDiagonalIndex(i);
        return idx >= 0 ? m_values.data() + idx : nullptr;
    }

    template <SpMatThreadSafety OtherThreadSafety>
    inline bool add(const SpMat<Scalar, Layout, Symmetry, OtherThreadSafety, Index>& other, Scalar w = 1.0)
    {
        assert(!m_filling && "add() must not be called between beginFill()/endFill()");
        assert(!other.isFilling() && "add() must not be called between beginFill()/endFill()");
        assert(other.rows() == m_rows && other.cols() == m_cols);

        const Index* const otherOuter = other.outerArray().data();
        const Index* const otherInner = other.innerArray().data();
        const Scalar* const otherValues = other.valuesArray().data();
        const Index* const innerPtr = m_inner.data();
        Scalar* const valuePtr = m_values.data();
        bool hasPending = false;
        for (Index o = 0; o < outerSize(); ++o)
        {
            Index p = m_outer[static_cast<size_t>(o)];
            const Index end = m_outer[static_cast<size_t>(o) + 1];
            for (Index q = otherOuter[o]; q < otherOuter[o + 1]; ++q)
            {
                const Index inner = otherInner[q];
                while (p < end && innerPtr[p] < inner)
                {
                    ++p;
                }
                if (p < end && innerPtr[p] == inner)
                {
                    valuePtr[p] += w * otherValues[q];
                }
                else
                {
                    m_buckets[static_cast<size_t>(o)].push_back(BucketEntry { inner, w * otherValues[q] });
                    hasPending = true;
                }
            }
        }
        if (hasPending)
        {
            mergeBuckets();
        }
        return hasPending;
    }

    inline MapType map()
    {
        assert(!m_filling && "map() called while filling");
        return MapType(m_rows, m_cols, nonZeros(), m_outer.data(), m_inner.data(), m_values.data());
    }

    inline ConstMapType map() const
    {
        assert(!m_filling && "map() called while filling");
        return ConstMapType(m_rows, m_cols, nonZeros(), m_outer.data(), m_inner.data(), m_values.data());
    }

    inline auto selfAdjointView()
    {
        static_assert(IsSymmetric, "selfAdjointView() requires a symmetric SpMat");
        return map().template selfadjointView<Symmetry == SpMatSymmetry::Upper ? Eigen::Upper : Eigen::Lower>();
    }

    inline auto selfAdjointView() const
    {
        static_assert(IsSymmetric, "selfAdjointView() requires a symmetric SpMat");
        return map().template selfadjointView<Symmetry == SpMatSymmetry::Upper ? Eigen::Upper : Eigen::Lower>();
    }

    template <class Derived>
    inline void setFromEigen(const Eigen::SparseMatrixBase<Derived>& expr)
    {
        assert(!m_filling && "setFromEigen() must not be called between beginFill()/endFill()");
        assert((!IsSymmetric || expr.rows() == expr.cols()) && "symmetric SpMat must be square");

        const auto assignCompressed = [this](const auto& src)
        {
            m_rows = static_cast<Index>(src.rows());
            m_cols = static_cast<Index>(src.cols());

            const Index* const outer = src.outerIndexPtr();
            const Index base = outer[0];
            const Index nnz = outer[outerSize()] - base;
            m_outer.resize(static_cast<size_t>(outerSize()) + 1);
            std::transform(outer, outer + outerSize() + 1, m_outer.begin(), [base](Index o) { return o - base; });
            m_inner.assign(src.innerIndexPtr() + base, src.innerIndexPtr() + base + nnz);
            m_values.assign(src.valuePtr() + base, src.valuePtr() + base + nnz);
            resetBucketStorage();
        };

        constexpr bool isCompatibleCompressed = !IsSymmetric && std::is_base_of_v<Eigen::SparseCompressedBase<Derived>, Derived> && std::is_same_v<typename Derived::Scalar, Scalar> &&
                                                std::is_same_v<typename Derived::StorageIndex, Index> && bool(Derived::IsRowMajor) == IsRowMajor;
        if constexpr (isCompatibleCompressed)
        {
            const Derived& src = expr.derived();
            if (src.isCompressed())
            {
                assignCompressed(src);
                return;
            }
        }

        EigenType src;
        if constexpr (IsSymmetric)
        {
            src = expr.derived().template triangularView<Symmetry == SpMatSymmetry::Upper ? Eigen::Upper : Eigen::Lower>();
        }
        else
        {
            src = expr;
        }
        src.makeCompressed();
        assignCompressed(src);
    }

private:
    static constexpr bool OuterStartsAtDiagonal = (Symmetry == SpMatSymmetry::Lower && IsColMajor) || (Symmetry == SpMatSymmetry::Upper && IsRowMajor);
    static constexpr bool OuterEndsAtDiagonal = (Symmetry == SpMatSymmetry::Upper && IsColMajor) || (Symmetry == SpMatSymmetry::Lower && IsRowMajor);

    inline Index findDiagonalIndex(Index i) const
    {
        assert(i >= 0 && i < m_rows && i < m_cols);
        const Index begin = m_outer[static_cast<size_t>(i)];
        const Index end = m_outer[static_cast<size_t>(i) + 1];
        if constexpr (OuterStartsAtDiagonal)
        {
            return begin < end && m_inner[static_cast<size_t>(begin)] == i ? begin : -1;
        }
        else if constexpr (OuterEndsAtDiagonal)
        {
            return begin < end && m_inner[static_cast<size_t>(end) - 1] == i ? end - 1 : -1;
        }
        else
        {
            return findIndexOuterInner(i, i);
        }
    }

    inline void resetBucketStorage()
    {
        m_buckets.clear();
        m_buckets.resize(static_cast<size_t>(outerSize()));
        if constexpr (IsThreadSafe)
        {
            m_outerMutexes = MutexArray(static_cast<size_t>(outerSize()));
        }
    }

    static constexpr bool HasAtomicAdd = std::is_same_v<Scalar, float> || std::is_same_v<Scalar, double> || (std::is_integral_v<Scalar> && !std::is_same_v<Scalar, bool>);

    inline std::unique_lock<std::mutex> makeOuterLock([[maybe_unused]] Index outer)
    {
        if constexpr (IsThreadSafe && !HasAtomicAdd)
        {
            return std::unique_lock<std::mutex>(m_outerMutexes[static_cast<size_t>(outer)]);
        }
        else
        {
            return std::unique_lock<std::mutex>();
        }
    }

    inline void ensureOuterLocked([[maybe_unused]] std::unique_lock<std::mutex>& lock, [[maybe_unused]] Index outer)
    {
        if constexpr (IsThreadSafe)
        {
            if (!lock.owns_lock())
            {
                lock = std::unique_lock<std::mutex>(m_outerMutexes[static_cast<size_t>(outer)]);
            }
        }
    }

    static inline void addValue(Scalar& target, Scalar value)
    {
        if constexpr (IsThreadSafe && HasAtomicAdd)
        {
            std::atomic_ref<Scalar>(target).fetch_add(value, std::memory_order_relaxed);
        }
        else
        {
            target += value;
        }
    }

    inline Index findIndex(Index row, Index col) const
    {
        assert(row >= 0 && row < m_rows && col >= 0 && col < m_cols);
        return findIndexOuterInner(IsRowMajor ? row : col, IsRowMajor ? col : row);
    }

    inline Index findIndexOuterInner(Index outer, Index inner) const
    {
        assert(outer >= 0 && outer < outerSize() && inner >= 0 && inner < innerSize());
        const Index begin = m_outer[static_cast<size_t>(outer)];
        const Index end = m_outer[static_cast<size_t>(outer) + 1];
        const Index* first = m_inner.data() + begin;
        const Index* last = m_inner.data() + end;
        const Index* it = std::lower_bound(first, last, inner);
        if (it == last || *it != inner)
        {
            return -1;
        }
        return static_cast<Index>(it - m_inner.data());
    }

    inline static void coalesce(Bucket& bucket)
    {
        if (bucket.size() < 2)
        {
            return;
        }
        std::sort(bucket.begin(), bucket.end(), [](const BucketEntry& a, const BucketEntry& b) { return a.m_inner < b.m_inner; });
        size_t w = 0;
        for (size_t r = 1; r < bucket.size(); ++r)
        {
            if (bucket[r].m_inner == bucket[w].m_inner)
            {
                bucket[w].m_value += bucket[r].m_value;
            }
            else
            {
                bucket[++w] = bucket[r];
            }
        }
        bucket.resize(w + 1);
    }

    inline void mergeBuckets()
    {
        const int nOuter = static_cast<int>(outerSize());

#pragma omp parallel for schedule(dynamic, 64)
        for (int o = 0; o < nOuter; ++o)
        {
            coalesce(m_buckets[static_cast<size_t>(o)]);
        }

        m_scratchOuter.resize(static_cast<size_t>(nOuter) + 1);
        m_scratchOuter[0] = Index(0);
        size_t total = 0;
        for (int o = 0; o < nOuter; ++o)
        {
            const size_t oldCount = static_cast<size_t>(m_outer[static_cast<size_t>(o) + 1] - m_outer[static_cast<size_t>(o)]);
            total += oldCount + m_buckets[static_cast<size_t>(o)].size();
            assert(total <= static_cast<size_t>(std::numeric_limits<Index>::max()) && "nnz overflows Index");
            m_scratchOuter[static_cast<size_t>(o) + 1] = static_cast<Index>(total);
        }
        m_scratchInner.resize(total);
        m_scratchValues.resize(total);

#pragma omp parallel for schedule(dynamic, 64)
        for (int o = 0; o < nOuter; ++o)
        {
            Bucket& bucket = m_buckets[static_cast<size_t>(o)];
            size_t p = static_cast<size_t>(m_outer[static_cast<size_t>(o)]);
            const size_t pEnd = static_cast<size_t>(m_outer[static_cast<size_t>(o) + 1]);
            size_t w = static_cast<size_t>(m_scratchOuter[static_cast<size_t>(o)]);

            if (bucket.empty())
            {
                std::copy(m_inner.begin() + p, m_inner.begin() + pEnd, m_scratchInner.begin() + w);
                std::copy(m_values.begin() + p, m_values.begin() + pEnd, m_scratchValues.begin() + w);
                continue;
            }

            for (const BucketEntry& e : bucket)
            {
                while (p < pEnd && m_inner[p] < e.m_inner)
                {
                    m_scratchInner[w] = m_inner[p];
                    m_scratchValues[w] = m_values[p];
                    ++p;
                    ++w;
                }
                assert((p == pEnd || m_inner[p] != e.m_inner) && "bucket entry already in pattern");
                m_scratchInner[w] = e.m_inner;
                m_scratchValues[w] = e.m_value;
                ++w;
            }
            std::copy(m_inner.begin() + p, m_inner.begin() + pEnd, m_scratchInner.begin() + w);
            std::copy(m_values.begin() + p, m_values.begin() + pEnd, m_scratchValues.begin() + w);
            bucket.clear();
        }

        m_outer.swap(m_scratchOuter);
        m_inner.swap(m_scratchInner);
        m_values.swap(m_scratchValues);
    }

    Index m_rows = 0;
    Index m_cols = 0;

    IndexArray m_outer;
    IndexArray m_inner;
    ScalarArray m_values;

    std::vector<Bucket> m_buckets;

    bool m_filling = false;
    IndexArray m_scratchOuter;
    IndexArray m_scratchInner;
    ScalarArray m_scratchValues;

    struct MutexArray
    {
        MutexArray() = default;
        explicit MutexArray(size_t n)
            : m_mutexes(n)
        {
        }
        MutexArray(const MutexArray& other)
            : m_mutexes(other.m_mutexes.size())
        {
        }
        MutexArray& operator=(const MutexArray& other)
        {
            if (this != &other)
            {
                m_mutexes = std::vector<std::mutex>(other.m_mutexes.size());
            }
            return *this;
        }
        MutexArray(MutexArray&&) noexcept = default;
        MutexArray& operator=(MutexArray&&) noexcept = default;

        inline size_t size() const
        {
            return m_mutexes.size();
        }
        inline std::mutex& operator[](size_t i)
        {
            return m_mutexes[i];
        }

        std::vector<std::mutex> m_mutexes;
    };

    MutexArray m_outerMutexes;
};

} // namespace spmat
