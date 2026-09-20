#pragma once

#include <cassert>
#include <cmath>

struct RobustSqrt
{
    virtual ~RobustSqrt() = default;

    virtual double evaluate(double x) const = 0;
    virtual double evaluateDeriv1(double x) const = 0;
    virtual double evaluateDeriv2(double x) const = 0;

    enum class Mode
    {
        Normal,
        Clamp
    };

    Mode m_mode = Mode::Normal;
};

struct SimpleRobustSqrt : public RobustSqrt
{
    explicit SimpleRobustSqrt(Mode mode, double eps = 1e-6)
        : m_eps(eps)
    {
        assert(m_eps > 0.0);
        m_mode = mode;
    }

    double evaluate(double x) const override
    {
        const double arg = x + m_eps;
        if (arg <= 0.0)
        {
            return m_mode == Mode::Clamp ? 0.0 : std::numeric_limits<double>::infinity();
        }
        return std::sqrt(arg) - std::sqrt(m_eps);
    }

    double evaluateDeriv1(double x) const override
    {
        const double arg = x + m_eps;
        if (arg <= 0.0)
        {
            return m_mode == Mode::Clamp ? 0.0 : std::numeric_limits<double>::infinity();
        }
        return 1.0 / (2.0 * std::sqrt(arg));
    }

    double evaluateDeriv2(double x) const override
    {
        const double arg = x + m_eps;
        if (arg <= 0.0)
        {
            return m_mode == Mode::Clamp ? 0.0 : std::numeric_limits<double>::infinity();
        }
        return -1.0 / (4.0 * std::pow(arg, 1.5));
    }

    double m_eps;
};