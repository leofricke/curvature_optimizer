#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <vector>

class SparseVector
{
public:
    struct Entry
    {
        int m_idx;
        double m_value;
    };

    void clear()
    {
        m_entries.clear();
    }

    const std::vector<Entry>& getEntries() const
    {
        return m_entries;
    }

    void addEntry(int idx, double value)
    {
        m_entries.emplace_back(idx, value);
    }

    void addSegment(int idx, const Eigen::Vector3d& segment)
    {
        m_entries.emplace_back(idx, segment[0]);
        m_entries.emplace_back(idx + 1, segment[1]);
        m_entries.emplace_back(idx + 2, segment[2]);
    }

    void sortEntries()
    {
        std::sort(m_entries.begin(), m_entries.end(), [](const Entry& a, const Entry& b) { return a.m_idx < b.m_idx; });
    }

private:
    std::vector<Entry> m_entries;
};