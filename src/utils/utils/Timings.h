#pragma once

#include "Timer.h"
#include <string>
#include <unordered_map>
#include <iomanip>

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
        if (!m_timers.contains(name))
        {
            m_timers[name] = Timer();
            m_measurements[name] = { 0.0 };
        }
        m_timers.at(name).start();
    }

    void stopTiming(const std::string& name)
    {
        m_timers.at(name).stop();

        m_measurements.at(name).push_back(m_timers.at(name).millis());
    }

    double lastMillis(const std::string& name) const
    {
        return m_measurements.at(name).back();
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
        double sum = 0.0;
        for (const double& m : m_measurements.at(name))
        {
            sum += m;
        }
        return sum / static_cast<double>(m_measurements.at(name).size());
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
