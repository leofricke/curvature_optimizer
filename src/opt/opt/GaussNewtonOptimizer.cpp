#include "GaussNewtonOptimizer.h"

#include "GaussNewtonEnergyFunction.h"
#include "LineSearch.h"
#include "LinearSolver.h"
#include "utils/SpMatUtil.h"
#include <utility>

static const double s_regularizationBeta = 1e-6;

GaussNewtonOptimizer::GaussNewtonOptimizer()
    : OptimizerBase(), m_dynamicInitialStepSize(false)
{
}

bool GaussNewtonOptimizer::getDynamicInitialStepSize() const
{
    return m_dynamicInitialStepSize;
}

void GaussNewtonOptimizer::setDynamicInitialStepSize(bool flag)
{
    m_dynamicInitialStepSize = flag;
}

double GaussNewtonOptimizer::run(std::shared_ptr<GaussNewtonEnergyFunction> energyFunc, Eigen::VectorXd& params)
{
    Eigen::VectorXd step;
    LinearSolver solver;

    step.resize(energyFunc->m_numVars);

    double prevEnergy = -1.0;
    double currentEnergy = energyFunc->evaluateResidualAndJacobian(params);
    onStart(params, currentEnergy);

    for (size_t i = 0; i < m_maxIterations; ++i)
    {
        m_statistics.m_timings.startTiming("Iteration Time");

        if (i > 0)
        {
            currentEnergy = energyFunc->evaluateResidualAndJacobian(params);
        }
        if (std::isnan(currentEnergy) || params.hasNaN())
        {
            m_statistics.m_status = Status::NanError;
            break;
        }

        const double gradNorm = energyFunc->m_system.m_gradient.norm();
        const double gradMaxNorm = energyFunc->m_system.m_gradient.lpNorm<Eigen::Infinity>();
        if (gradMaxNorm < m_gradientTolerance)
        {
            m_statistics.m_status = Status::GradientNormBelowTolerance;
            break;
        }

        double regularizer = 0.0;
        GaussNewtonSystem& system = energyFunc->m_system;
        const double minDiag = minDiagonal(system.m_hessian);
        if (minDiag <= 1e-12)
        {
            regularizer = -minDiag + s_regularizationBeta;
        }
        regularize(system.m_hessian, regularizer);

        const bool patternDirty = system.m_hessianPatternDirty || i == 0;
        if (!solver.compute(std::as_const(system.m_hessian).map(), patternDirty))
        {
            m_statistics.m_status = Status::SolverFailed;
            break;
        }
        if (!solver.solve(-energyFunc->m_system.m_gradient, step))
        {
            m_statistics.m_status = Status::SolverFailed;
            break;
        }
        if (step.dot(energyFunc->m_system.m_gradient) >= 0.0)
        {
            m_statistics.m_status = Status::NoDescentDirection;
            break;
        }

        const double stepNorm = step.norm();
        const double stepMaxNorm = step.lpNorm<Eigen::Infinity>();
        if (stepMaxNorm < m_stepTolerance)
        {
            m_statistics.m_status = Status::StepNormBelowTolerance;
            break;
        }

        double stepSize = 1.0;
        if (m_dynamicInitialStepSize)
        {
            if (i == 0)
            {
                stepSize = std::clamp(1.01 / gradMaxNorm, 1e-8, 1.0);
            }
            else
            {
                stepSize = std::clamp(2.02 * (currentEnergy - prevEnergy) / step.dot(energyFunc->m_system.m_gradient), 1e-8, 1.0);
            }
        }
        prevEnergy = currentEnergy;
        if (!lineSearchArmijo(params, currentEnergy, stepSize, step, *energyFunc))
        {
            m_statistics.m_status = Status::NoImprovement;
            break;
        }

        m_statistics.m_timings.stopTiming("Iteration Time");

        onStep(params, i, currentEnergy, stepSize, stepNorm, gradNorm);
        m_statistics.m_numIterations++;

        if (prevEnergy - currentEnergy < m_energyReductionTolerance)
        {
            m_statistics.m_status = Status::EnergyReductionBelowTolerance;
            break;
        }
    }
    m_statistics.m_timings.stopTiming("Iteration Time");

    if (m_statistics.m_status == Status::Running)
    {
        m_statistics.m_status = Status::MaxIterationsExhausted;
    }
    onStop(params, currentEnergy);

    return currentEnergy;
}
