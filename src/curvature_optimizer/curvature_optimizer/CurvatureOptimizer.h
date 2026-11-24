#pragma once

#include "IVertexCurvatureOperator.h"
#include "quantities/DualAreaDerivatives.h"
#include "TargetData.h"
#include "utils/Barycentric.h"
#include "utils/Timings.h"

#include <pmp/surface_mesh.h>
#include <Eigen/Core>
#include <map>

struct CurvatureOptimizerWorkspace;
struct ModularGaussNewtonEnergyFunction;

/**
* Main algorithm class for curvature optimization.
* It allows to set target principal curvatures, target boundary curvatures, and position constraints for mesh vertices.
* The mesh geometry is optimized to achieve these targets while minimizing additional energy terms to preserve area and conformality.
* Supports both AutoForm::Mesh and pmp::SurfaceMesh as input and output.
*/
class CurvatureOptimizer
{
public:
    /**
    * TargetCurvature: optimize for target principal curvatures at vertices.
    * WillmoreEnergy: optimize for Willmore energy of the mesh.
    * CurvatureVariation: optimize for variation of curvature over the mesh (does not need target curvatures).
    */
    enum class OptimizerMode
    {
        TargetCurvature = 0,
        WillmoreEnergy = 1,
        CurvatureVariation = 2
    };

    /**
    * Status of one optimization step, including the current energy and the energies of each individual blocks.
    */
    struct Progress
    {
        double m_energy;
        std::map<std::string, double> m_blockEnergies;
    };

    using MeshUpdateCallback = std::function<void(const Progress& progress)>;

    explicit CurvatureOptimizer(std::shared_ptr<pmp::SurfaceMesh> mesh);
    virtual ~CurvatureOptimizer();

    void clearTargets(bool clearMetric = true);
    void setTargetMetric(const pmp::SurfaceMesh& referenceMesh);

    void addTargetVertexCurvature(pmp::Vertex v, double kmax, double kmin);
    void addTargetBoundaryCurvature(pmp::Vertex v, double kb);
    void addTargetPosition(const BarycentricPoint& bary, const Eigen::Vector3d& position);

    const std::vector<TargetVertexCurvature>& getTargetVertexCurvatures() const;
    const std::vector<TargetBoundaryCurvature>& getTargetBoundaryCurvatures() const;
    const std::vector<TargetPosition>& getTargetPositions() const;

    double getCurvatureEnergyWeight() const;
    double getBoundaryCurvatureEnergyWeight() const;
    double getAreaEnergyWeight() const;
    double getConformalEnergyWeight() const;
    double getPositionEnergyWeight() const;
    void setCurvatureEnergyWeight(double w);
    void setBoundaryCurvatureEnergyWeight(double w);
    void setAreaEnergyWeight(double w);
    void setConformalEnergyWeight(double w);
    void setPositionEnergyWeight(double w);

    OptimizerMode getOptimizerMode() const;
    void setOptimizerMode(OptimizerMode mode);
    bool getDynamicInitialStepSize() const;
    void setDynamicInitialStepSize(bool flag);
    size_t getMaxIterations() const;
    void setMaxIterations(size_t maxIterations);
    double getPNormFairing() const;
    void setPNormFairing(double pNorm);

    bool getUseRobustDihedralAngles() const;
    void setUseRobustDihedralAngles(bool flag);
    DualAreaMode getDualAreaMode() const;
    void setDualAreaMode(DualAreaMode mode);
    VertexCurvatureOperator getVertexCurvatureOperator() const;
    void setVertexCurvatureOperator(VertexCurvatureOperator op);

    /**
    * Enables verbose output during optimization.
    */
    void setVerbose(bool verbose);
    /**
    * Callback invoked after each optimization step.
    */
    void setMeshUpdateCallback(MeshUpdateCallback callback);

private:
    void initTargetMetric();
    std::shared_ptr<ModularGaussNewtonEnergyFunction> createEnergyFunction() const;

public:
    /**
    * Must be called before invoking compute().
    * After initMesh() was called, this class may also be used to compute current curvatures using for example getCurrentVertexCurvature(...) .
    */
    void initMesh();
    /**
    * Perform the curvature optimization. Locked vertices are kept fixed using hard constraints.
    */
    Timings compute(const std::vector<pmp::Vertex>& lockedVertices = {});

    const PrincipalCurvatures& getCurrentVertexCurvature(pmp::Vertex v) const;
    double getCurrentBoundaryCurvature(pmp::Vertex v) const;

    double getCurvatureRelativeError() const;
    size_t getNumberOfExtremeDihedralAngles() const;

private:
    std::shared_ptr<pmp::SurfaceMesh> m_mesh;

    OptimizerMode m_optimizerMode;
    bool m_dynamicInitialStepSize;
    size_t m_maxIterations;
    double m_pNormFairing;
    double m_curvatureEnergyWeight;
    double m_boundaryCurvatureEnergyWeight;
    double m_areaEnergyWeight;
    double m_conformalEnergyWeight;
    double m_positionEnergyWeight;

    std::vector<TargetVertexCurvature> m_targetVertexCurvatures;
    std::vector<TargetBoundaryCurvature> m_targetBoundaryCurvatures;
    std::vector<TargetMetric> m_targetMetrics;
    std::vector<TargetPosition> m_targetPositions;

    std::unique_ptr<CurvatureOptimizerWorkspace> m_ws;

    double m_initialCurvatureError;
    double m_currentCurvatureError;

    bool m_verbose;
    MeshUpdateCallback m_meshUpdateCallback;
};