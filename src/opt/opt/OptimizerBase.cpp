#include "OptimizerBase.h"

#include <iostream>
#include <iomanip>

OptimizerBase::OptimizerBase()
    : m_maxIterations(500), m_gradientTolerance(1e-8), m_stepTolerance(1e-8), m_energyReductionTolerance(1e-8), m_constantMatrixPattern(false), m_verbose(false)
{
}

void OptimizerBase::setMaxIterations(size_t maxIterations)
{
    m_maxIterations = maxIterations;
}

void OptimizerBase::setGradientTolerance(double tolerance)
{
    m_gradientTolerance = tolerance;
}

void OptimizerBase::setStepTolerance(double tolerance)
{
    m_stepTolerance = tolerance;
}

void OptimizerBase::setEnergyReductionTolerance(double tolerance)
{
    m_energyReductionTolerance = tolerance;
}

void OptimizerBase::setConstantMatrixPattern(bool flag)
{
    m_constantMatrixPattern = flag;
}

void OptimizerBase::setVerbose(bool flag)
{
    m_verbose = flag;
}

void OptimizerBase::setUpdateCallback(UpdateCallback callback)
{
    m_updateCallback = callback;
}

const OptimizerBase::Statistics& OptimizerBase::getStatistics() const
{
    return m_statistics;
}

void OptimizerBase::onStart(const Eigen::VectorXd& params, double currentEnergy)
{
    if (m_updateCallback)
    {
        m_updateCallback(params, currentEnergy);
    }

    m_statistics = Statistics();
    m_statistics.m_status = Status::Running;
    m_statistics.m_history.clear();
    m_statistics.m_history.reserve(m_maxIterations);
    m_statistics.m_history.emplace_back(currentEnergy);
    m_statistics.m_timings.startTiming("Total Time");

    if (m_verbose)
    {
        printHeader(currentEnergy);
    }
}

void OptimizerBase::onStep(const Eigen::VectorXd& params, size_t iteration, double currentEnergy, double stepSize, double stepNorm, double gradNorm)
{
    if (m_updateCallback)
    {
        m_updateCallback(params, currentEnergy);
    }
    m_statistics.m_history.emplace_back(currentEnergy);

    if (m_verbose)
    {
        printIteration(iteration, currentEnergy, stepSize, stepNorm, gradNorm);
    }
}

void OptimizerBase::onStop(const Eigen::VectorXd& /*params*/, double /*currentEnergy*/)
{
    m_statistics.m_timings.stopTiming("Total Time");

    if (m_verbose)
    {
        std::cout << m_statistics << std::endl;
    }
}

void OptimizerBase::printHeader(double currentEnergy) const
{
    const int precision = 5;
    const int spacing = 10 + precision;
    std::cout << std::setw(4) << "Iter";
    std::cout << std::setw(10) << "Time";
    std::cout << std::setw(spacing) << "Energy";
    std::cout << std::setw(spacing) << "StepSize";
    std::cout << std::setw(spacing) << "StepNorm";
    std::cout << std::setw(spacing) << "GradNorm";
    std::cout << std::endl;

    printIteration(-1, currentEnergy, NAN, NAN, NAN);
}

void OptimizerBase::printIteration(int iteration, double currentEnergy, double stepSize, double stepNorm, double gradNorm) const
{
    constexpr int precision = 5;
    constexpr int spacing = 10 + precision;
    if (iteration >= 0)
        std::cout << std::setw(4) << iteration;
    else
        std::cout << std::setw(4) << "";

    if (iteration >= 0)
        std::cout << std::setw(10) << m_statistics.m_timings.lastTime("Iteration Time");
    else
        std::cout << std::setw(10) << "";

    if (!std::isnan(currentEnergy))
        std::cout << std::setw(spacing) << std::scientific << std::setprecision(precision) << currentEnergy;
    else
        std::cout << std::setw(spacing) << "";

    if (!std::isnan(stepSize))
        std::cout << std::setw(spacing) << std::scientific << std::setprecision(precision) << stepSize;
    else
        std::cout << std::setw(spacing) << "";

    if (!std::isnan(stepNorm))
        std::cout << std::setw(spacing) << std::scientific << std::setprecision(precision) << stepNorm;
    else
        std::cout << std::setw(spacing) << "";

    if (!std::isnan(gradNorm))
        std::cout << std::setw(spacing) << std::scientific << std::setprecision(precision) << gradNorm;
    else
        std::cout << std::setw(spacing) << "";

    std::cout << std::endl << std::fixed;
}

std::ostream& operator<<(std::ostream& os, const OptimizerBase::Statistics& stats)
{
    std::string status = "";
    switch (stats.m_status)
    {
    case OptimizerBase::Status::NotStarted:
        status = "NotStarted";
        break;
    case OptimizerBase::Status::Running:
        status = "Running";
        break;
    case OptimizerBase::Status::MaxIterationsExhausted:
        status = "MaxIterationsExhausted";
        break;
    case OptimizerBase::Status::NoImprovement:
        status = "NoImprovement";
        break;
    case OptimizerBase::Status::GradientNormBelowTolerance:
        status = "GradientNormBelowTolerance";
        break;
    case OptimizerBase::Status::StepNormBelowTolerance:
        status = "StepNormBelowTolerance";
        break;
    case OptimizerBase::Status::EnergyReductionBelowTolerance:
        status = "EnergyReductionBelowTolerance";
        break;
    case OptimizerBase::Status::SolverFailed:
        status = "SolverFailed";
        break;
    case OptimizerBase::Status::NoDescentDirection:
        status = "NoDescentDirection";
        break;
    case OptimizerBase::Status::NanError:
        status = "NanError";
        break;
    }
    const double finalEnergy = stats.m_history.empty() ? NAN : stats.m_history.back().m_energy;
    os << "Result:    FinalEnergy: " << finalEnergy << ", Status: " << status << ", NumIterations: " << stats.m_numIterations << std::endl
       << "Timing:    TotalTime: " << stats.m_timings.accumulatedTime("Total Time") << ", MeanIterationTime: " << stats.m_timings.meanTime("Iteration Time");

    return os;
}
