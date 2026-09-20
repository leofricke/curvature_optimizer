#pragma once

#include <Eigen/Core>

struct EnergyFunction;

/**
* Performs a backtracking line-search to find a step size that satisfies the Armijo condition.
*/
bool lineSearchArmijo(Eigen::VectorXd& params, double& currentEnergy, double& stepSize, const Eigen::VectorXd& stepDir, EnergyFunction& energyFunc);