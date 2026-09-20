#pragma once

#include "IVertexCurvatureOperator.h"
#include "CurvatureTensorVertexCurvature.h"
#include "MeanGaussianVertexCurvature.h"
#include "quantities/DualAreaDerivatives.h"
#include "opt/EnergyFunction.h"
#include "opt/MeshDofs.h"

#include <pmp/surface_mesh.h>
#include <variant>

/**
* Acts as the main data store during optimization. Stores the current mesh geometry, basic mesh properties, vertex curvatures, and their derivatives.
*/
struct CurvatureOptimizerWorkspace
{
    using VertexCurvatureVariantType = std::variant<CurvatureTensorVertexCurvature, MeanGaussianVertexCurvature>;

    explicit CurvatureOptimizerWorkspace(pmp::SurfaceMesh& mesh);
    virtual ~CurvatureOptimizerWorkspace();

    void initMeshData();
    void reset();

    void updateVertexCurvatureOperator();

    void updateGeometry(EnergyFunction::DerivativeMode mode = EnergyFunction::DerivativeMode::NoDerivs);
    Eigen::Vector3d position(pmp::Vertex v) const;

    const IVertexCurvatureOperator& getCurrentVertexCurvature(pmp::Vertex v) const;
    IVertexCurvatureOperator& getCurrentVertexCurvature(pmp::Vertex v);

    pmp::SurfaceMesh& m_mesh;

    bool m_initialized;
    MeshDofs m_dofs;

    double m_currentMeshArea;
    pmp::FaceProperty<double> m_currentFaceAreas;
    pmp::FaceProperty<Eigen::Vector3d> m_currentFaceNormals;
    pmp::HalfedgeProperty<double> m_currentInteriorAngles;
    pmp::EdgeProperty<double> m_currentDihedralAngles;
    pmp::VertexProperty<double> m_currentVertexDualAreas;
    pmp::VertexProperty<VertexCurvatureVariantType> m_currentVertexCurvatures;

    pmp::FaceProperty<LocalGradient<Eigen::RowVector3d, 4>> m_currentFaceAreaGradients;
    pmp::HalfedgeProperty<LocalGradient<Eigen::RowVector3d, 3>> m_currentInteriorAngleGradients;
    pmp::EdgeProperty<LocalGradient<Eigen::RowVector3d, 4>> m_currentDihedralAngleGradients;
    pmp::VertexProperty<LocalGradient<Eigen::RowVector3d, 8>> m_currentVertexDualAreaGradients;

    pmp::VertexProperty<LocalGradient<Eigen::RowVector3d, 8>> m_currentMaxCurvatureGradients;
    pmp::VertexProperty<LocalGradient<Eigen::RowVector3d, 8>> m_currentMinCurvatureGradients;

    bool m_robustDihedralAngles;
    DualAreaMode m_dualAreaMode;
    VertexCurvatureOperator m_vertexCurvatureOperator;
};