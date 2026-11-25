#pragma once

#include "Timer.h"
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

class Timings
{
public:
    Timings() = default;

    void reset()
    {
        m_timers.clear();
        m_measurements.clear();
    }

    void startTiming(const std::string& name)
    {
        if (auto timer = m_timers.find(name); timer != m_timers.end())
        {
            timer->second.start();
        }
        else
        {
            m_timers[name] = Timer(true);
            m_measurements[name] = {};
        }
    }

    void stopTiming(const std::string& name)
    {
        const auto it = m_timers.find(name);
        if (it == m_timers.end() || !it->second.running())
        {
            return;
        }
        it->second.stop();
        m_measurements.at(name).push_back(it->second.millis());
    }

    double lastMillis(const std::string& name) const
    {
        const std::vector<double>& measurements = m_measurements.at(name);
        return measurements.empty() ? 0.0 : measurements.back();
    }

    double accumulatedMillis(const std::string& name) const
    {
        double sum = 0.0;
        for (const double& m : m_measurements.at(name))
        {
            sum += m;
        }
        return sum;
    }

    double meanMillis(const std::string& name) const
    {
        const std::vector<double>& measurements = m_measurements.at(name);
        if (measurements.empty())
        {
            return 0.0;
        }
        return accumulatedMillis(name) / static_cast<double>(measurements.size());
    }

    static std::string formatMillis(double millis, double precision = 1)
    {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision) << millis << " ms";
        return oss.str();
    }

    std::string lastTime(const std::string& name, double precision = 1) const
    {
        return formatMillis(lastMillis(name), precision);
    }

    std::string accumulatedTime(const std::string& name, double precision = 1) const
    {
        return formatMillis(accumulatedMillis(name), precision);
    }

    std::string meanTime(const std::string& name, double precision = 1) const
    {
        return formatMillis(meanMillis(name), precision);
    }

private:
    std::unordered_map<std::string, Timer> m_timers;
    std::unordered_map<std::string, std::vector<double>> m_measurements;
};
