#pragma once

#include <Eigen/Core>

struct PrincipalCurvatures
{
    PrincipalCurvatures()
        : m_kmax(0.0), m_kmin(0.0)
    {
        m_maxDir.setZero();
        m_minDir.setZero();
    }

    PrincipalCurvatures(double kmax, double kmin, const Eigen::Vector3d& maxDir, const Eigen::Vector3d& minDir)
        : m_kmax(kmax), m_kmin(kmin), m_maxDir(maxDir), m_minDir(minDir)
    {
    }

    double maxAbsCurvature() const
    {
        return std::max(std::abs(m_kmax), std::abs(m_kmin));
    }

    double meanCurvature() const
    {
        return 0.5 * (m_kmax + m_kmin);
    }

    double gaussianCurvature() const
    {
        return m_kmax * m_kmin;
    }

    double m_kmax;
    double m_kmin;
    Eigen::Vector3d m_maxDir;
    Eigen::Vector3d m_minDir;
};