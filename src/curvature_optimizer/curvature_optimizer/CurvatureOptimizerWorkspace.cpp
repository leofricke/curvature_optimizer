#include "CurvatureOptimizerWorkspace.h"

#include "quantities/DihedralAngleDerivatives.h"
#include "quantities/InteriorAngleDerivatives.h"
#include "quantities/FaceAreaDerivatives.h"

CurvatureOptimizerWorkspace::CurvatureOptimizerWorkspace(pmp::SurfaceMesh& mesh)
    : m_mesh(mesh),
      m_initialized(false),
      m_dofs(m_mesh),
      m_currentMeshArea(0.0),
      m_robustDihedralAngles(true),
      m_dualAreaMode(DualAreaMode::VoronoiMixed),
      m_vertexCurvatureOperator(VertexCurvatureOperator::MeanGaussian)
{
}

CurvatureOptimizerWorkspace::~CurvatureOptimizerWorkspace()
{
    reset();
}

void CurvatureOptimizerWorkspace::initMeshData()
{
    m_currentFaceAreas = m_mesh.face_property<double>("opt:currentFaceAreas", 0.0);
    m_currentFaceNormals = m_mesh.face_property<Eigen::Vector3d>("opt:currentFaceNormals", Eigen::Vector3d::Zero());
    m_currentInteriorAngles = m_mesh.halfedge_property<double>("opt:currentInteriorAngles", 0.0);
    m_currentDihedralAngles = m_mesh.edge_property<double>("opt:currentDihedralAngles", 0.0);
    m_currentVertexDualAreas = m_mesh.vertex_property<double>("opt:currentVertexDualAreas", 0.0);
    m_currentVertexCurvatures = m_mesh.vertex_property<VertexCurvatureVariantType>("opt:currentVertexCurvatures");

    m_currentFaceAreaGradients = m_mesh.face_property<LocalGradient<Eigen::RowVector3d, 4>>("opt:currentFaceAreaGradients");
    m_currentInteriorAngleGradients = m_mesh.halfedge_property<LocalGradient<Eigen::RowVector3d, 3>>("opt:currentInteriorAngleGradients");
    m_currentDihedralAngleGradients = m_mesh.edge_property<LocalGradient<Eigen::RowVector3d, 4>>("opt:currentDihedralAngleGradients");
    m_currentVertexDualAreaGradients = m_mesh.vertex_property<LocalGradient<Eigen::RowVector3d, 8>>("opt:currentVertexDualAreaGradients");

    m_currentMaxCurvatureGradients = m_mesh.vertex_property<LocalGradient<Eigen::RowVector3d, 8>>("opt:currentMaxCurvatureGradients");
    m_currentMinCurvatureGradients = m_mesh.vertex_property<LocalGradient<Eigen::RowVector3d, 8>>("opt:currentMinCurvatureGradients");

    updateVertexCurvatureOperator();
    updateGeometry();

    m_initialized = true;
}

void CurvatureOptimizerWorkspace::reset()
{
    m_initialized = false;

    m_dofs.reset();

    m_currentFaceAreas = m_mesh.get_face_property<double>("opt:currentFaceAreas");
    m_currentFaceNormals = m_mesh.get_face_property<Eigen::Vector3d>("opt:currentFaceNormals");
    m_currentInteriorAngles = m_mesh.get_halfedge_property<double>("opt:currentInteriorAngles");
    m_currentDihedralAngles = m_mesh.get_edge_property<double>("opt:currentDihedralAngles");
    m_currentVertexDualAreas = m_mesh.get_vertex_property<double>("opt:currentVertexDualAreas");
    m_currentVertexCurvatures = m_mesh.get_vertex_property<VertexCurvatureVariantType>("opt:currentVertexCurvatures");

    m_currentFaceAreaGradients = m_mesh.get_face_property<LocalGradient<Eigen::RowVector3d, 4>>("opt:currentFaceAreaGradients");
    m_currentInteriorAngleGradients = m_mesh.get_halfedge_property<LocalGradient<Eigen::RowVector3d, 3>>("opt:currentInteriorAngleGradients");
    m_currentDihedralAngleGradients = m_mesh.get_edge_property<LocalGradient<Eigen::RowVector3d, 4>>("opt:currentDihedralAngleGradients");
    m_currentVertexDualAreaGradients = m_mesh.get_vertex_property<LocalGradient<Eigen::RowVector3d, 8>>("opt:currentVertexDualAreaGradients");

    m_currentMaxCurvatureGradients = m_mesh.get_vertex_property<LocalGradient<Eigen::RowVector3d, 8>>("opt:currentMaxCurvatureGradients");
    m_currentMinCurvatureGradients = m_mesh.get_vertex_property<LocalGradient<Eigen::RowVector3d, 8>>("opt:currentMinCurvatureGradients");

    if (m_currentFaceAreas)
    {
        m_mesh.remove_face_property<double>(m_currentFaceAreas);
    }
    if (m_currentFaceNormals)
    {
        m_mesh.remove_face_property<Eigen::Vector3d>(m_currentFaceNormals);
    }
    if (m_currentInteriorAngles)
    {
        m_mesh.remove_halfedge_property<double>(m_currentInteriorAngles);
    }
    if (m_currentDihedralAngles)
    {
        m_mesh.remove_edge_property<double>(m_currentDihedralAngles);
    }
    if (m_currentVertexDualAreas)
    {
        m_mesh.remove_vertex_property<double>(m_currentVertexDualAreas);
    }
    if (m_currentVertexCurvatures)
    {
        m_mesh.remove_vertex_property<VertexCurvatureVariantType>(m_currentVertexCurvatures);
    }
    if (m_currentFaceAreaGradients)
    {
        m_mesh.remove_face_property<LocalGradient<Eigen::RowVector3d, 4>>(m_currentFaceAreaGradients);
    }
    if (m_currentInteriorAngleGradients)
    {
        m_mesh.remove_halfedge_property<LocalGradient<Eigen::RowVector3d, 3>>(m_currentInteriorAngleGradients);
    }
    if (m_currentDihedralAngleGradients)
    {
        m_mesh.remove_edge_property<LocalGradient<Eigen::RowVector3d, 4>>(m_currentDihedralAngleGradients);
    }
    if (m_currentVertexDualAreaGradients)
    {
        m_mesh.remove_vertex_property<LocalGradient<Eigen::RowVector3d, 8>>(m_currentVertexDualAreaGradients);
    }
    if (m_currentMaxCurvatureGradients)
    {
        m_mesh.remove_vertex_property<LocalGradient<Eigen::RowVector3d, 8>>(m_currentMaxCurvatureGradients);
    }
    if (m_currentMinCurvatureGradients)
    {
        m_mesh.remove_vertex_property<LocalGradient<Eigen::RowVector3d, 8>>(m_currentMinCurvatureGradients);
    }
}

void CurvatureOptimizerWorkspace::updateVertexCurvatureOperator()
{
    if (m_vertexCurvatureOperator == VertexCurvatureOperator::CurvatureTensor)
    {
        for (pmp::Vertex v : m_mesh.vertices())
        {
            m_currentVertexCurvatures[v] = CurvatureTensorVertexCurvature();
        }
    }
    else if (m_vertexCurvatureOperator == VertexCurvatureOperator::MeanGaussian)
    {
        for (pmp::Vertex v : m_mesh.vertices())
        {
            m_currentVertexCurvatures[v] = MeanGaussianVertexCurvature();
        }
    }
    else
    {
        throw std::runtime_error("Invalid curvature operator");
    }
    updateGeometry();
}

void CurvatureOptimizerWorkspace::updateGeometry(EnergyFunction::DerivativeMode mode)
{
    // update areas and normals
    m_currentMeshArea = 0.0;
    for (pmp::Face f : m_mesh.faces())
    {
        m_currentFaceNormals[f] = computeFaceNormal(*this, f);

        if (mode == EnergyFunction::DerivativeMode::NoDerivs)
        {
            m_currentFaceAreas[f] = computeFaceArea(*this, f);
        }
        else if (mode == EnergyFunction::DerivativeMode::FirstDerivs)
        {
            m_currentFaceAreaGradients[f].clear();
            m_currentFaceAreas[f] = computeFaceAreaGradient(*this, f, m_currentFaceAreaGradients[f]);
        }

        m_currentMeshArea += m_currentFaceAreas[f];
    }

    // update interior angles
    for (pmp::Halfedge h : m_mesh.halfedges())
    {
        if (mode == EnergyFunction::DerivativeMode::NoDerivs)
        {
            m_currentInteriorAngles[h] = computeInteriorAngle(*this, h);
        }
        else if (mode == EnergyFunction::DerivativeMode::FirstDerivs)
        {
            m_currentInteriorAngleGradients[h].clear();
            m_currentInteriorAngles[h] = computeInteriorAngleGradient(*this, h, m_currentInteriorAngleGradients[h]);
        }
    }

    // update dihedral angles
    for (pmp::Edge e : m_mesh.edges())
    {
        if (mode == EnergyFunction::DerivativeMode::NoDerivs)
        {
            m_currentDihedralAngles[e] = computeDihedralAngle(*this, e, m_robustDihedralAngles);
        }
        else if (mode == EnergyFunction::DerivativeMode::FirstDerivs)
        {
            m_currentDihedralAngleGradients[e].clear();
            m_currentDihedralAngles[e] = computeDihedralAngleGradient(*this, e, m_currentDihedralAngleGradients[e], m_robustDihedralAngles);
        }
    }

    // update vertex dual areas
    for (pmp::Vertex v : m_mesh.vertices())
    {
        if (mode == EnergyFunction::DerivativeMode::NoDerivs)
        {
            m_currentVertexDualAreas[v] = computeDualArea(*this, v, m_dualAreaMode);
        }
        else if (mode == EnergyFunction::DerivativeMode::FirstDerivs)
        {
            m_currentVertexDualAreaGradients[v].clear();
            m_currentVertexDualAreas[v] = computeDualAreaGradient(*this, v, m_currentVertexDualAreaGradients[v], m_dualAreaMode);
        }
    }

    // update curvatures
    for (pmp::Vertex v : m_mesh.vertices())
    {
        IVertexCurvatureOperator& currentCurv = getCurrentVertexCurvature(v);

        if (mode == EnergyFunction::DerivativeMode::NoDerivs)
        {
            currentCurv.computePrincipalCurvatures(*this, v);
        }
        else if (mode == EnergyFunction::DerivativeMode::FirstDerivs)
        {
            m_currentMaxCurvatureGradients[v].clear();
            m_currentMinCurvatureGradients[v].clear();
            currentCurv.computePrincipalCurvaturesGradient(*this, v, m_currentMaxCurvatureGradients[v], m_currentMinCurvatureGradients[v]);
        }
    }
}

Eigen::Vector3d CurvatureOptimizerWorkspace::position(pmp::Vertex v) const
{
    return m_dofs.position(v);
}

const IVertexCurvatureOperator& CurvatureOptimizerWorkspace::getCurrentVertexCurvature(pmp::Vertex v) const
{
    if (m_vertexCurvatureOperator == VertexCurvatureOperator::CurvatureTensor)
    {
        return std::get<CurvatureTensorVertexCurvature>(m_currentVertexCurvatures[v]);
    }
    else if (m_vertexCurvatureOperator == VertexCurvatureOperator::MeanGaussian)
    {
        return std::get<MeanGaussianVertexCurvature>(m_currentVertexCurvatures[v]);
    }
    else
    {
        throw std::runtime_error("Invalid curvature operator");
    }
}

IVertexCurvatureOperator& CurvatureOptimizerWorkspace::getCurrentVertexCurvature(pmp::Vertex v)
{
    if (m_vertexCurvatureOperator == VertexCurvatureOperator::CurvatureTensor)
    {
        return std::get<CurvatureTensorVertexCurvature>(m_currentVertexCurvatures[v]);
    }
    else if (m_vertexCurvatureOperator == VertexCurvatureOperator::MeanGaussian)
    {
        return std::get<MeanGaussianVertexCurvature>(m_currentVertexCurvatures[v]);
    }
    else
    {
        throw std::runtime_error("Invalid curvature operator");
    }
}
