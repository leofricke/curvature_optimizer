#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <type_traits>
#include <vector>

namespace spmat
{

template <typename Scalar, typename Index = int>
class SpVec
{
    static_assert(std::is_integral_v<Index> && std::is_signed_v<Index>, "IndexType must be a signed integral type");

public:
    struct Entry
    {
        Index m_idx;
        Scalar m_value;
    };

    SpVec()
        : SpVec(0)
    {
    }

    explicit SpVec(Index size)
        : m_size(size)
    {
        assert(size >= 0);
    }

    SpVec(Index size, const std::vector<Entry>& entries)
    {
        setFromEntries(size, entries);
    }

    inline void resize(Index size)
    {
        assert(!m_filling && "resize() must not be called between beginFill()/endFill()");
        assert(size >= 0);
        if (size == m_size)
        {
            return;
        }
        m_size = size;
        m_entries.clear();
    }

    inline void reserve(Index nnz)
    {
        assert(!m_filling && "reserve() must not be called between beginFill()/endFill()");
        assert(nnz >= 0);
        m_entries.reserve(static_cast<std::size_t>(nnz));
    }

    inline void clearPattern()
    {
        assert(!m_filling && "clearPattern() must not be called between beginFill()/endFill()");
        m_entries.clear();
    }

    inline bool isFilling() const
    {
        return m_filling;
    }

    inline void beginFill()
    {
        assert(!m_filling && "beginFill() called twice without endFill()");
        setZero();
        m_fillStart = static_cast<Index>(m_entries.size());
        m_filling = true;
    }

    inline bool endFill()
    {
        assert(m_filling && "endFill() without a matching beginFill()");
        m_filling = false;

        if (static_cast<Index>(m_entries.size()) == m_fillStart)
        {
            return false;
        }
        sortAndMerge(m_fillStart);
        return true;
    }

    inline void setZero()
    {
        for (Entry& entry : m_entries)
        {
            entry.m_value = Scalar(0.0);
        }
    }

    inline Index size() const
    {
        return m_size;
    }

    inline Index nonZeros() const
    {
        return static_cast<Index>(m_entries.size());
    }

    inline const std::vector<Entry>& entries() const
    {
        return m_entries;
    }

    inline void setFromEntries(Index size, const std::vector<Entry>& entries)
    {
        assert(!m_filling && "setFromEntries() must not be called between beginFill()/endFill()");
        assert(size >= 0);
        assert(std::all_of(entries.begin(), entries.end(), [size](const Entry& e) { return e.m_idx >= 0 && e.m_idx < size; }) && "setFromEntries() entry index out of range");
        m_size = size;
        m_entries = entries;
        if (m_entries.empty())
        {
            return;
        }
        sortAndMerge(0);
    }

    inline void insert(Index idx, Scalar value)
    {
        assert(idx >= 0 && idx < m_size);
        if (m_filling)
        {
            m_entries.push_back(Entry { idx, value });
            return;
        }
        const auto it = std::lower_bound(m_entries.begin(), m_entries.end(), idx, [](const Entry& e, Index i) { return e.m_idx < i; });

        if (it != m_entries.end() && it->m_idx == idx)
        {
            it->m_value += value;
            return;
        }
        m_entries.insert(it, Entry { idx, value });
    }

    template <typename Derived>
    inline void insertSegment(Index idx, const Derived& segment)
    {
        const Index segmentSize = static_cast<Index>(segment.size());
        assert(idx >= 0 && idx + segmentSize <= m_size);
        if (segmentSize == 0)
        {
            return;
        }
        if (m_filling)
        {
            for (Index k = 0; k < segmentSize; ++k)
            {
                m_entries.push_back(Entry { idx + k, segment[k] });
            }
            return;
        }
        const auto first = std::lower_bound(m_entries.begin(), m_entries.end(), idx, [](const Entry& e, Index i) { return e.m_idx < i; });
        const auto firstPos = static_cast<Index>(first - m_entries.begin());

        const auto last = std::lower_bound(first, m_entries.end(), idx + segmentSize, [](const Entry& e, Index i) { return e.m_idx < i; });

        const Index existing = static_cast<Index>(last - first);
        const Index missing = segmentSize - existing;
        if (missing > 0)
        {
            m_entries.insert(last, static_cast<std::size_t>(missing), Entry { 0, Scalar(0) });
        }

        Index read = firstPos + existing - 1;
        for (Index k = segmentSize - 1; k >= 0; --k)
        {
            Scalar value = segment[k];
            if (read >= firstPos && m_entries[read].m_idx == idx + k)
            {
                value += m_entries[read].m_value;
                --read;
            }
            Entry& entry = m_entries[firstPos + k];
            entry.m_idx = idx + k;
            entry.m_value = value;
        }
    }

    inline Scalar coeff(Index idx) const
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        assert(idx >= 0 && idx < m_size);
        const Index p = findIndex(idx);
        return p >= 0 ? m_entries[p].m_value : Scalar(0);
    }

    inline Scalar& coeffRef(Index idx)
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        assert(idx >= 0 && idx < m_size);
        const Index p = findIndex(idx);
        assert(p >= 0 && "coeffRef() on a coefficient that is not in the pattern");
        return m_entries[p].m_value;
    }

    inline Scalar* coeffPtr(Index idx)
    {
        assert(!m_filling && "Coefficient access is not allowed while filling");
        assert(idx >= 0 && idx < m_size);
        const Index p = findIndex(idx);
        return p >= 0 ? &m_entries[p].m_value : nullptr;
    }

private:
    inline void sortAndMerge(Index sortedPrefix)
    {
        const auto byIdx = [](const Entry& a, const Entry& b)
        {
            return a.m_idx < b.m_idx;
        };
        const auto mid = m_entries.begin() + sortedPrefix;
        std::sort(mid, m_entries.end(), byIdx);
        std::inplace_merge(m_entries.begin(), mid, m_entries.end(), byIdx);

        auto write = m_entries.begin();
        for (auto read = std::next(m_entries.begin()); read != m_entries.end(); ++read)
        {
            if (read->m_idx == write->m_idx)
            {
                write->m_value += read->m_value;
            }
            else
            {
                ++write;
                *write = *read;
            }
        }
        m_entries.erase(std::next(write), m_entries.end());
    }

    inline Index findIndex(Index idx) const
    {
        const auto begin = m_entries.begin();
        const auto end = m_entries.end();
        const auto it = std::lower_bound(begin, end, idx, [](const Entry& e, Index i) { return e.m_idx < i; });
        if (it == end || it->m_idx != idx)
        {
            return -1;
        }
        return static_cast<Index>(it - begin);
    }

    Index m_size;
    std::vector<Entry> m_entries;
    bool m_filling = false;
    Index m_fillStart = 0;
};

} // namespace spmat
