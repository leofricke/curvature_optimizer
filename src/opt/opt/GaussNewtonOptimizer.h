#pragma once

#include "OptimizerBase.h"
#include <memory>

struct GaussNewtonEnergyFunction;

/**
* Gauss-Newton optimization algorithm.
* Requires a energy function of type GaussNewtonEnergyFunction that can evaluate the energy and build the system.
* The optimizer iteratively solves linear systems to obtain a step direction and then performs a line-search to find steps that decrease the energy.
*/
class GaussNewtonOptimizer : public OptimizerBase
{
public:
    GaussNewtonOptimizer();

    bool getDynamicInitialStepSize() const;
    void setDynamicInitialStepSize(bool flag);

    double run(std::shared_ptr<GaussNewtonEnergyFunction> energyFunc, Eigen::VectorXd& params);

private:
    bool m_dynamicInitialStepSize;
};