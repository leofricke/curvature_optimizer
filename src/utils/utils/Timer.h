#pragma once

#include <chrono>
#include <iostream>
#include <stdexcept>

class Timer
{
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;

public:
    Timer(bool autostart = false)
        : m_running(false)
    {
        if (autostart)
        {
            start();
        }
    }

    bool running() const
    {
        return m_running;
    }

    void start()
    {
        if (m_running)
        {
            throw std::runtime_error("Timer already running");
        }
        m_startTime = Clock::now();
        m_running = true;
    }

    void stop()
    {
        if (!m_running)
        {
            throw std::runtime_error("Timer not running");
        }
        m_stopTime = Clock::now();
        m_running = false;
    }

    double millis() const
    {
        if (m_running)
        {
            throw std::runtime_error("Timer still running");
        }
        const Duration timeSpan = std::chrono::duration_cast<Duration>(m_stopTime - m_startTime);
        return 1000.0 * timeSpan.count();
    }

private:
    TimePoint m_startTime;
    TimePoint m_stopTime;
    bool m_running;
};

inline std::ostream& operator<<(std::ostream& os, const Timer& timer)
{
    const auto p = os.precision();
    os.precision(2);
    os << timer.millis() << " ms";
    os.precision(p);
    return os;
}