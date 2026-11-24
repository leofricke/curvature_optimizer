#include <gtest/gtest.h>

#include "curvature_optimizer/VirtualSubdivision.h"

TEST(VirtualSubdivisionTest, BasicTestSingleTriangle)
{
    pmp::SurfaceMesh mesh;
    const pmp::Vertex v0 = mesh.add_vertex(pmp::Point(0.0, 0.0, 0.0));
    const pmp::Vertex v1 = mesh.add_vertex(pmp::Point(1.0, 0.0, 0.0));
    const pmp::Vertex v2 = mesh.add_vertex(pmp::Point(1.0, 2.0, 0.0));
    const pmp::Face f = mesh.add_triangle(v0, v1, v2);

    const pmp::Halfedge h = mesh.halfedge(f);
    const Eigen::Vector3d p0 = mesh.position(mesh.from_vertex(h));
    const Eigen::Vector3d p1 = mesh.position(mesh.to_vertex(h));
    const Eigen::Vector3d p2 = mesh.position(mesh.to_vertex(mesh.next_halfedge(h)));

    VirtualSubdivision subdiv;
    subdiv.updateCenterWeights(mesh, f);

    InlineVector<VirtualTriangle, 4> triangles;
    subdiv.initVirtualTriangles(mesh, f, triangles);
    const VirtualTriangle& tri = triangles[0];

    const Eigen::Vector3d incenter = tri.m_incenter.interpolate(p0, p1, p2);
    const Eigen::Vector3d q01 = tri.m_q01.interpolate(p0, p1, p2);
    const Eigen::Vector3d q12 = tri.m_q12.interpolate(p0, p1, p2);
    const Eigen::Vector3d q20 = tri.m_q20.interpolate(p0, p1, p2);

    EXPECT_TRUE(tri.m_incenter.isInside());

    // test qxx lie on edges
    EXPECT_NEAR((q01 - p0).normalized().dot((p1 - p0).normalized()), 1.0, 1e-6);
    EXPECT_NEAR((q12 - p1).normalized().dot((p2 - p1).normalized()), 1.0, 1e-6);
    EXPECT_NEAR((q20 - p2).normalized().dot((p0 - p2).normalized()), 1.0, 1e-6);

    // test incircle radius
    const double d0 = (q01 - incenter).squaredNorm();
    const double d1 = (q12 - incenter).squaredNorm();
    const double d2 = (q20 - incenter).squaredNorm();
    EXPECT_NEAR(d0, d1, 1e-6);
    EXPECT_NEAR(d1, d2, 1e-6);
    EXPECT_NEAR(d2, d0, 1e-6);

    // test bisected (outer) angles
    const double alpha00 = std::acos(((incenter - p0).normalized().dot((p1 - p0).normalized())));
    const double alpha01 = std::acos(((incenter - p0).normalized().dot((p2 - p0).normalized())));
    const double alpha10 = std::acos(((incenter - p1).normalized().dot((p0 - p1).normalized())));
    const double alpha11 = std::acos(((incenter - p1).normalized().dot((p2 - p1).normalized())));
    const double alpha20 = std::acos(((incenter - p2).normalized().dot((p0 - p2).normalized())));
    const double alpha21 = std::acos(((incenter - p2).normalized().dot((p1 - p2).normalized())));
    EXPECT_NEAR(alpha00, alpha01, 1e-6);
    EXPECT_NEAR(alpha10, alpha11, 1e-6);
    EXPECT_NEAR(alpha20, alpha21, 1e-6);

    // test bisected (inner) angles
    const double beta00 = std::acos(((q01 - incenter).normalized().dot((p0 - incenter).normalized())));
    const double beta01 = std::acos(((q20 - incenter).normalized().dot((p0 - incenter).normalized())));
    const double beta10 = std::acos(((q01 - incenter).normalized().dot((p1 - incenter).normalized())));
    const double beta11 = std::acos(((q12 - incenter).normalized().dot((p1 - incenter).normalized())));
    const double beta20 = std::acos(((q12 - incenter).normalized().dot((p2 - incenter).normalized())));
    const double beta21 = std::acos(((q20 - incenter).normalized().dot((p2 - incenter).normalized())));
    EXPECT_NEAR(beta00, beta01, 1e-6);
    EXPECT_NEAR(beta10, beta11, 1e-6);
    EXPECT_NEAR(beta20, beta21, 1e-6);
}