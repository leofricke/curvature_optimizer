#pragma once

#include <pmp/surface_mesh.h>
#include <Eigen/Core>
#include <vector>

/**
* Maps vertex positions from a mesh to a flat DoF vector.
*/
class MeshDofs
{
public:
    MeshDofs(pmp::SurfaceMesh& mesh)
        : m_mesh(mesh)
    {
        reset();
    }

    void reset()
    {
        m_vertToDofIdx.clear();
        m_dofs.resize(0);
    }

    void initMeshDofs(const std::vector<pmp::Vertex>& lockedVertices)
    {
        m_vertToDofIdx.assign(m_mesh.vertices_size(), 0);
        for (pmp::Vertex locked : lockedVertices)
        {
            m_vertToDofIdx[locked.idx()] = -1;
        }

        size_t idx = 0;
        for (pmp::Vertex v : m_mesh.vertices())
        {
            if (m_vertToDofIdx[v.idx()] >= 0)
            {
                m_vertToDofIdx[v.idx()] = idx++;
            }
        }

        m_dofs.resize(idx * 3);
        for (pmp::Vertex v : m_mesh.vertices())
        {
            const int dofIdx = m_vertToDofIdx[v.idx()];
            if (dofIdx >= 0)
            {
                m_dofs.segment<3>(dofIdx * 3) = static_cast<Eigen::Vector3d>(m_mesh.position(v));
            }
        }
    }

    int dofIdx(pmp::Vertex v) const
    {
        return m_vertToDofIdx[v.idx()];
    }

    size_t numDofs() const
    {
        return m_dofs.size();
    }

    const Eigen::VectorXd& getDofs() const
    {
        return m_dofs;
    }

    Eigen::VectorXd& getDofs()
    {
        return m_dofs;
    }

    Eigen::Vector3d position(pmp::Vertex v) const
    {
        if (m_dofs.size() > 0)
        {
            const int dofIdx = m_vertToDofIdx[v.idx()];
            if (dofIdx >= 0)
            {
                return m_dofs.segment<3>(dofIdx * 3);
            }
        }

        return static_cast<Eigen::Vector3d>(m_mesh.position(v));
    }

    void updateMesh()
    {
        for (pmp::Vertex v : m_mesh.vertices())
        {
            const int dofIdx = m_vertToDofIdx[v.idx()];
            if (dofIdx >= 0)
            {
                m_mesh.position(v) = m_dofs.segment<3>(dofIdx * 3);
            }
        }
    }

private:
    pmp::SurfaceMesh& m_mesh;
    std::vector<int> m_vertToDofIdx;
    Eigen::VectorXd m_dofs;
};