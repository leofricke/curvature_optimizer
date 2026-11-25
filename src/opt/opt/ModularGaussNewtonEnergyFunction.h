#pragma once

#include "GaussNewtonEnergyFunction.h"
#include "GaussNewtonSystem.h"
#include <memory>

/**
* A Gauss-Newton energy function for use in a larger modular energy function.
* Acts similar to a regular GaussNewtonEnergyFunction, but does not own a residual vector or any parts of the Gauss-Newton system.
* In the evaluateResidual and evaluateResidualAndJacobian methods, the block may only access the residual vector and Jacobian rows from m_offset to m_offset + m_numResiduals.
* The m_offset value is set externally by the ModularGaussNewtonEnergyFunction.
*/
struct ModularGaussNewtonEnergyBlock
{
    ModularGaussNewtonEnergyBlock(size_t numResiduals, size_t numVars)
        : m_numResiduals(numResiduals), m_numVars(numVars), m_offset(0), m_name("name not set")
    {
    }

    virtual ~ModularGaussNewtonEnergyBlock() = default;

    virtual void evaluateResidual(const Eigen::VectorXd& params, GaussNewtonSystem& system) = 0;
    virtual void evaluateResidualAndJacobian(const Eigen::VectorXd& params, GaussNewtonSystem& system) = 0;

    /**
    * Called by the ModularGaussNewtonEnergyFunction to initialize weights of this block into the global weight vector.
    */
    virtual void initWeights(double energyWeight, Eigen::VectorXd& weights) = 0;

    size_t m_numResiduals;
    size_t m_numVars;
    size_t m_offset;
    std::string m_name;
};

/**
* Represents a Gauss-Newton energy function which is a sum of multiple individual energies.
* Each individual energy must be provided as a ModularGaussNewtonEnergyBlock, which has methods to evaluate its residuals and Jacobians.
* In contrast to EnergyFunction, the blocks do not own a residual vector or any parts of the Gauss-Newton system.
* Once the ModularGaussNewtonEnergyFunction with a given set of blocks, they are stacked into a large residual vector and the m_offset of each block is set accordingly.
*/
struct ModularGaussNewtonEnergyFunction : public GaussNewtonEnergyFunction
{
    ModularGaussNewtonEnergyFunction(const std::vector<std::shared_ptr<ModularGaussNewtonEnergyBlock>>& blocks)
        : GaussNewtonEnergyFunction(0, 0), m_blocks(blocks)
    {
        m_numVars = m_blocks[0]->m_numVars;

        size_t offset = 0;
        for (auto& block : m_blocks)
        {
            if (block->m_numVars != m_numVars)
            {
                throw std::runtime_error("Number of variables must be equal for all blocks");
            }

            block->m_offset = offset;
            offset += block->m_numResiduals;
            m_numResiduals += block->m_numResiduals;
        }

        m_system.init(m_numResiduals, m_numVars);
    }

    double evaluateResidualImpl(const Eigen::VectorXd& params) override
    {
        for (auto& block : m_blocks)
        {
            block->evaluateResidual(params, m_system);
        }
        return m_system.evaluateEnergy();
    }

    double evaluateResidualAndJacobianImpl(const Eigen::VectorXd& params) override
    {
        m_system.beginFill();
        for (auto& block : m_blocks)
        {
            block->evaluateResidualAndJacobian(params, m_system);
        }
        m_system.endFill();
        return m_system.evaluateEnergy();
    }

    /**
    * Calls initWeights on each block. Accepts a block weight vector to weigh the blocks differently.
    */
    void initWeights(const std::vector<double>& blockWeights)
    {
        if (m_blocks.size() != blockWeights.size())
        {
            throw std::runtime_error("Weights array has wrong size");
        }

        for (size_t i = 0; i < m_blocks.size(); ++i)
        {
            m_blocks[i]->initWeights(blockWeights[i], m_system.m_weights);
        }
    }

    /**
    * Evaluate the energy value of the block at given index.
    */
    double getBlockEnergy(size_t idx) const
    {
        const auto& block = m_blocks[idx];
        // clang-format off
        return 0.5 * m_system.m_residual.segment(block->m_offset, block->m_numResiduals).dot(
            m_system.m_weights.segment(block->m_offset, block->m_numResiduals).asDiagonal() *
            m_system.m_residual.segment(block->m_offset, block->m_numResiduals));
        // clang-format on
    }

    /**
    * Evaluates the energy values of all blocks.
    */
    std::map<std::string, double> getBlockEnergies() const
    {
        std::map<std::string, double> blockEnergies;

        for (size_t i = 0; i < m_blocks.size(); ++i)
        {
            blockEnergies.emplace(m_blocks[i]->m_name, getBlockEnergy(i));
        }

        return blockEnergies;
    }

    std::vector<std::shared_ptr<ModularGaussNewtonEnergyBlock>> m_blocks;
};
