#include "CurvatureOptimizer.h"

#include "CurvatureOptimizerWorkspace.h"
#include "PrincipalCurvatures.h"
#include "energies/AreaEnergy.h"
#include "energies/BoundaryCurvatureEnergy.h"
#include "energies/ConformalEnergy.h"
#include "energies/CurvatureEnergy.h"
#include "energies/CurvatureVariationEnergy.h"
#include "energies/PolygonPlanarityEnergy.h"
#include "energies/PositionEnergy.h"
#include "energies/RegularizationEnergy.h"
#include "energies/WillmoreEnergy.h"
#include "quantities/BoundaryCurvatureDerivatives.h"
#include "opt/GaussNewtonOptimizer.h"
#include "opt/GaussNewtonEnergyFunction.h"
#include <numbers>
#include <stdexcept>

CurvatureOptimizer::CurvatureOptimizer(std::shared_ptr<pmp::SurfaceMesh> mesh)
    : m_mesh(mesh),
      m_optimizerMode(OptimizerMode::TargetCurvature),
      m_dynamicInitialStepSize(true),
      m_maxIterations(100),
      m_pNormFairing(2.0),
      m_curvatureEnergyWeight(1.0),
      m_boundaryCurvatureEnergyWeight(1.0),
      m_areaEnergyWeight(1.0),
      m_conformalEnergyWeight(1.0),
      m_positionEnergyWeight(1.0),
      m_initialCurvatureError(0.0),
      m_currentCurvatureError(0.0),
      m_verbose(false)
{
    m_ws = std::make_unique<CurvatureOptimizerWorkspace>(*m_mesh);
}

CurvatureOptimizer::~CurvatureOptimizer()
{
}

void CurvatureOptimizer::clearTargets(bool clearMetric)
{
    m_targetVertexCurvatures.clear();
    m_targetBoundaryCurvatures.clear();
    m_targetPositions.clear();
    if (clearMetric)
    {
        m_targetMetrics.clear();
    }
}

void CurvatureOptimizer::setTargetMetric(const pmp::SurfaceMesh& referenceMesh)
{
    if (referenceMesh.vertices_size() != m_mesh->vertices_size() || referenceMesh.halfedges_size() != m_mesh->halfedges_size() || referenceMesh.faces_size() != m_mesh->faces_size())
    {
        throw std::invalid_argument("setTargetMetric: reference mesh must have the same topology");
    }

    m_targetMetrics.clear();
    for (pmp::Face f : referenceMesh.faces())
    {
        m_targetMetrics.push_back(createTargetMetric(referenceMesh, f));
    }
}

void CurvatureOptimizer::addTargetVertexCurvature(pmp::Vertex v, double kmax, double kmin)
{
    m_targetVertexCurvatures.emplace_back(v, kmax, kmin);
}

void CurvatureOptimizer::addTargetBoundaryCurvature(pmp::Vertex v, double kb)
{
    if (!m_mesh->is_boundary(v))
    {
        throw std::runtime_error("Not a boundary vertex");
    }
    m_targetBoundaryCurvatures.emplace_back(v, kb);
}

void CurvatureOptimizer::addTargetPosition(const BarycentricPoint& bary, const Eigen::Vector3d& position)
{
    m_targetPositions.emplace_back(bary, position);
}

const std::vector<TargetVertexCurvature>& CurvatureOptimizer::getTargetVertexCurvatures() const
{
    return m_targetVertexCurvatures;
}

const std::vector<TargetBoundaryCurvature>& CurvatureOptimizer::getTargetBoundaryCurvatures() const
{
    return m_targetBoundaryCurvatures;
}

const std::vector<TargetPosition>& CurvatureOptimizer::getTargetPositions() const
{
    return m_targetPositions;
}

double CurvatureOptimizer::getCurvatureEnergyWeight() const
{
    return m_curvatureEnergyWeight;
}

double CurvatureOptimizer::getBoundaryCurvatureEnergyWeight() const
{
    return m_boundaryCurvatureEnergyWeight;
}

double CurvatureOptimizer::getAreaEnergyWeight() const
{
    return m_areaEnergyWeight;
}

double CurvatureOptimizer::getConformalEnergyWeight() const
{
    return m_conformalEnergyWeight;
}

double CurvatureOptimizer::getPositionEnergyWeight() const
{
    return m_positionEnergyWeight;
}

void CurvatureOptimizer::setCurvatureEnergyWeight(double w)
{
    m_curvatureEnergyWeight = w;
}

void CurvatureOptimizer::setBoundaryCurvatureEnergyWeight(double w)
{
    m_boundaryCurvatureEnergyWeight = w;
}

void CurvatureOptimizer::setAreaEnergyWeight(double w)
{
    m_areaEnergyWeight = w;
}

void CurvatureOptimizer::setConformalEnergyWeight(double w)
{
    m_conformalEnergyWeight = w;
}

void CurvatureOptimizer::setPositionEnergyWeight(double w)
{
    m_positionEnergyWeight = w;
}

CurvatureOptimizer::OptimizerMode CurvatureOptimizer::getOptimizerMode() const
{
    return m_optimizerMode;
}

void CurvatureOptimizer::setOptimizerMode(OptimizerMode mode)
{
    m_optimizerMode = mode;
}

bool CurvatureOptimizer::getDynamicInitialStepSize() const
{
    return m_dynamicInitialStepSize;
}

void CurvatureOptimizer::setDynamicInitialStepSize(bool flag)
{
    m_dynamicInitialStepSize = flag;
}

size_t CurvatureOptimizer::getMaxIterations() const
{
    return m_maxIterations;
}

void CurvatureOptimizer::setMaxIterations(size_t maxIterations)
{
    m_maxIterations = maxIterations;
}

double CurvatureOptimizer::getPNormFairing() const
{
    return m_pNormFairing;
}

void CurvatureOptimizer::setPNormFairing(double pNorm)
{
    m_pNormFairing = pNorm;
}

bool CurvatureOptimizer::getUseRobustDihedralAngles() const
{
    return m_ws->m_robustDihedralAngles;
}

void CurvatureOptimizer::setUseRobustDihedralAngles(bool flag)
{
    m_ws->m_robustDihedralAngles = flag;
}

DualAreaMode CurvatureOptimizer::getDualAreaMode() const
{
    return m_ws->m_dualAreaMode;
}

void CurvatureOptimizer::setDualAreaMode(DualAreaMode mode)
{
    m_ws->m_dualAreaMode = mode;
}

VertexCurvatureOperator CurvatureOptimizer::getVertexCurvatureOperator() const
{
    return m_ws->m_vertexCurvatureOperator;
}

void CurvatureOptimizer::setVertexCurvatureOperator(VertexCurvatureOperator op)
{
    m_ws->m_vertexCurvatureOperator = op;
    if (m_ws->m_initialized)
    {
        m_ws->updateVertexCurvatureOperator();
    }
}

void CurvatureOptimizer::setVerbose(bool verbose)
{
    m_verbose = verbose;
}

void CurvatureOptimizer::setMeshUpdateCallback(MeshUpdateCallback callback)
{
    m_meshUpdateCallback = callback;
}

void CurvatureOptimizer::initTargetMetric()
{
    if (!m_targetMetrics.empty())
    {
        return;
    }

    for (pmp::Face f : m_mesh->faces())
    {
        m_targetMetrics.push_back(createTargetMetric(*m_mesh, f));
    }
}

std::shared_ptr<GaussNewtonEnergyFunction> CurvatureOptimizer::createEnergyFunction() const
{
    std::vector<std::shared_ptr<GaussNewtonEnergyModule>> modules;
    std::vector<double> moduleWeights;
    if (m_optimizerMode == OptimizerMode::TargetCurvature)
    {
        modules.push_back(std::make_shared<CurvatureEnergy>(*m_ws, m_targetVertexCurvatures));
        moduleWeights.push_back(m_curvatureEnergyWeight);
        if (!m_targetBoundaryCurvatures.empty())
        {
            modules.push_back(std::make_shared<BoundaryCurvatureEnergy>(*m_ws, m_targetBoundaryCurvatures));
            moduleWeights.push_back(m_boundaryCurvatureEnergyWeight);
        }
    }
    else if (m_optimizerMode == OptimizerMode::WillmoreEnergy)
    {
        modules.push_back(std::make_shared<WillmoreEnergy>(*m_ws));
        moduleWeights.push_back(m_curvatureEnergyWeight);
    }
    else if (m_optimizerMode == OptimizerMode::CurvatureVariation)
    {
        modules.push_back(std::make_shared<CurvatureVariationEnergy>(*m_ws, m_pNormFairing));
        moduleWeights.push_back(m_curvatureEnergyWeight);
    }
    modules.push_back(std::make_shared<AreaEnergy>(*m_ws, m_targetMetrics));
    moduleWeights.push_back(m_areaEnergyWeight);
    modules.push_back(std::make_shared<ConformalEnergy>(*m_ws, m_targetMetrics));
    moduleWeights.push_back(m_conformalEnergyWeight);
    if (!m_mesh->is_triangle_mesh())
    {
        modules.push_back(std::make_shared<PolygonPlanarityEnergy>(*m_ws));
        moduleWeights.push_back(1.0);
    }
    if (!m_targetPositions.empty())
    {
        modules.push_back(std::make_shared<PositionEnergy>(*m_ws, m_targetPositions));
        moduleWeights.push_back(m_positionEnergyWeight);
    }
    modules.push_back(std::make_shared<RegularizationEnergy>(*m_ws));
    moduleWeights.push_back(1e-3);

    std::shared_ptr<GaussNewtonEnergyFunction> energyFunc = std::make_shared<GaussNewtonEnergyFunction>(modules, moduleWeights);

    return energyFunc;
}

void CurvatureOptimizer::initMesh()
{
    m_ws->initMeshData();

    m_initialCurvatureError = 0.0;
    m_currentCurvatureError = 0.0;
}

Timings CurvatureOptimizer::compute(const std::vector<pmp::Vertex>& lockedVertices)
{
    if (!m_ws->m_initialized)
    {
        throw std::runtime_error("Call initMesh() first");
    }

    m_ws->updateGeometry();
    initTargetMetric();

    m_ws->m_dofs.initMeshDofs(lockedVertices);

    std::shared_ptr<GaussNewtonEnergyFunction> energyFunc = createEnergyFunction();

    m_initialCurvatureError = 0.0;
    m_currentCurvatureError = 0.0;

    GaussNewtonOptimizer optimizer;
    optimizer.setMaxIterations(m_maxIterations);
    optimizer.setDynamicInitialStepSize(m_dynamicInitialStepSize);
    optimizer.setVerbose(m_verbose);

    energyFunc->addBeforeEvaluateCallback([&](const Eigen::VectorXd&, EnergyFunction::DerivativeMode mode) { m_ws->updateGeometry(mode); });

    optimizer.setUpdateCallback(
        [&](const Eigen::VectorXd&, double currentEnergy)
        {
            m_ws->m_dofs.updateMesh();

            if (m_initialCurvatureError <= 0.0)
            {
                m_initialCurvatureError = m_currentCurvatureError = energyFunc->getModuleEnergy(0);
            }
            else
            {
                m_currentCurvatureError = energyFunc->getModuleEnergy(0);
            }

            if (m_meshUpdateCallback)
            {
                const Progress progress = { currentEnergy, energyFunc->getModuleEnergies() };
                m_meshUpdateCallback(progress);
            }
        });

    optimizer.run(energyFunc, m_ws->m_dofs.getDofs());

    m_ws->m_dofs.updateMesh();
    m_ws->updateGeometry();
    m_ws->m_dofs.reset();

    return optimizer.getStatistics().m_timings;
}

const PrincipalCurvatures& CurvatureOptimizer::getCurrentVertexCurvature(pmp::Vertex v) const
{
    return m_ws->getCurrentVertexCurvature(v).m_curv;
}

double CurvatureOptimizer::getCurrentBoundaryCurvature(pmp::Vertex v) const
{
    if (!m_ws->m_mesh.is_boundary(v))
    {
        throw std::runtime_error("Not a boundary vertex");
    }

    return computeBoundaryCurvature(*m_ws, v);
}

double CurvatureOptimizer::getCurvatureRelativeError() const
{
    if (m_initialCurvatureError <= 0.0)
    {
        return 0.0;
    }

    return m_currentCurvatureError / m_initialCurvatureError;
}

size_t CurvatureOptimizer::getNumberOfExtremeDihedralAngles() const
{
    size_t n = 0;
    for (pmp::Edge e : m_mesh->edges())
    {
        const double beta = m_ws->m_currentDihedralAngles[e];
        if (std::abs(beta) > 0.75 * std::numbers::pi)
        {
            n++;
        }
    }

    return n;
}