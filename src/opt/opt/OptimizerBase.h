#pragma once

#include "utils/Timings.h"
#include <Eigen/Core>

/**
* Base class for non-linear optimization algorithms.
* Provides common functionality for running an optimization, convergence conditions, and statistics tracking.
*/
class OptimizerBase
{
public:
    using UpdateCallback = std::function<void(const Eigen::VectorXd& params, double currentEnergy)>;

    enum class Status
    {
        NotStarted,
        Running,
        MaxIterationsExhausted,
        NoImprovement,
        GradientNormBelowTolerance,
        StepNormBelowTolerance,
        EnergyReductionBelowTolerance,
        SolverFailed,
        NoDescentDirection,
        NanError
    };

    struct StepStatistics
    {
        double m_energy;
    };
    struct Statistics
    {
        size_t m_numIterations = 0;
        Status m_status = Status::NotStarted;
        Timings m_timings;
        std::vector<StepStatistics> m_history = {};
    };

    OptimizerBase();

    void setMaxIterations(size_t maxIterations);
    void setGradientTolerance(double tolerance);
    void setStepTolerance(double tolerance);
    void setEnergyReductionTolerance(double tolerance);
    void setConstantMatrixPattern(bool flag);
    void setVerbose(bool flag);

    void setUpdateCallback(UpdateCallback callback);
    const Statistics& getStatistics() const;

protected:
    void onStart(const Eigen::VectorXd& params, double currentEnergy);
    void onStep(const Eigen::VectorXd& params, size_t iteration, double currentEnergy, double stepSize, double stepNorm, double gradNorm);
    void onStop(const Eigen::VectorXd& params, double currentEnergy);

private:
    void printHeader(double currentEnergy) const;
    void printIteration(int iteration, double currentEnergy, double stepSize, double stepNorm, double gradNorm) const;

protected:
    size_t m_maxIterations;
    double m_gradientTolerance;
    double m_stepTolerance;
    double m_energyReductionTolerance;
    bool m_constantMatrixPattern;
    bool m_verbose;

    UpdateCallback m_updateCallback;
    Statistics m_statistics;
};

std::ostream& operator<<(std::ostream& os, const OptimizerBase::Statistics& stats);
