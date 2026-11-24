#include "LineSearch.h"

#include "EnergyFunction.h"
#include <Eigen/Core>

static const double s_lineSearchMinAlpha = 1e-12;
static const double s_armijoCoeff = 1e-4; // armijo sufficient decrease coefficient

struct LineSearchSample
{
    double m_alpha;
    double m_energy;
    double m_gradient;
    bool m_hasGradient;
};

bool lineSearchArmijo(Eigen::VectorXd& params, double& currentEnergy, double& stepSize, const Eigen::VectorXd& stepDir, EnergyFunction& energyFunc)
{
    const Eigen::VectorXd origin = params;

    LineSearchSample initial;
    initial.m_alpha = 0.0;
    initial.m_energy = currentEnergy;
    initial.m_gradient = stepDir.dot(energyFunc.gradient());
    initial.m_hasGradient = true;

    LineSearchSample current;
    current.m_alpha = stepSize > 0.0 ? stepSize : 1.0;
    current.m_hasGradient = false;
    while (true)
    {
        params = origin + current.m_alpha * stepDir;
        current.m_energy = energyFunc.evaluateEnergy(params);
        current.m_hasGradient = false;

        const bool valid = !std::isnan(current.m_energy);
        if (valid && current.m_energy <= initial.m_energy + s_armijoCoeff * current.m_alpha * initial.m_gradient)
        {
            currentEnergy = current.m_energy;
            stepSize = current.m_alpha;

            return true;
        }

        current.m_alpha *= 0.5;

        if (current.m_alpha < s_lineSearchMinAlpha)
        {
            break;
        }
    }
    params = origin;

    return false;
}
