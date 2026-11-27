#include <gtest/gtest.h>

#include "curvature_optimizer/CurvatureOptimizer.h"

#include "TestFilesUtil.h"
#include <pmp/io/io.h>

class CurvatureOptimizerSphereTest : public ::testing::Test
{
public:
    CurvatureOptimizerSphereTest()
    {
        mesh = std::make_shared<pmp::SurfaceMesh>();
        read(*mesh, getTestFile("sphere.off"));
    }

    std::shared_ptr<pmp::SurfaceMesh> mesh;
};

class CurvatureOptimizerHemisphereTest : public ::testing::Test
{
public:
    CurvatureOptimizerHemisphereTest()
    {
        mesh = std::make_shared<pmp::SurfaceMesh>();
        read(*mesh, getTestFile("hemisphere.off"));
    }

    std::shared_ptr<pmp::SurfaceMesh> mesh;
};

class CurvatureOptimizerCylinderTest : public ::testing::Test
{
public:
    CurvatureOptimizerCylinderTest()
    {
        mesh = std::make_shared<pmp::SurfaceMesh>();
        read(*mesh, getTestFile("cylinder.off"));
    }

    std::shared_ptr<pmp::SurfaceMesh> mesh;
};

TEST_F(CurvatureOptimizerSphereTest, ScaleSphere)
{
    CurvatureOptimizer opt(mesh);
    opt.setAreaEnergyWeight(0.0);
    opt.setConformalEnergyWeight(1.0);
    opt.initMesh();
    for (pmp::Vertex v : mesh->vertices())
    {
        opt.addTargetVertexCurvature(v, 1.0, 1.0);
    }
    opt.compute();

    for (pmp::Vertex v : mesh->vertices())
    {
        const auto& curv = opt.getCurrentVertexCurvature(v);
        EXPECT_NEAR(curv.m_kmax, 1.0, 0.1);
        EXPECT_NEAR(curv.m_kmin, 1.0, 0.1);
    }
}

TEST_F(CurvatureOptimizerHemisphereTest, Plane)
{
    CurvatureOptimizer opt(mesh);
    opt.initMesh();
    for (pmp::Vertex v : mesh->vertices())
    {
        opt.addTargetVertexCurvature(v, 0.0, 0.0);
    }
    opt.compute();

    for (pmp::Vertex v : mesh->vertices())
    {
        const auto& curv = opt.getCurrentVertexCurvature(v);
        EXPECT_NEAR(curv.m_kmax, 0.0, 0.1);
        EXPECT_NEAR(curv.m_kmin, 0.0, 0.1);
    }
}

TEST_F(CurvatureOptimizerHemisphereTest, Invert)
{
    CurvatureOptimizer opt(mesh);
    opt.initMesh();
    for (pmp::Vertex v : mesh->vertices())
    {
        opt.addTargetVertexCurvature(v, -1.0, -1.0);
    }
    opt.compute();

    for (pmp::Vertex v : mesh->vertices())
    {
        if (!mesh->is_boundary(v))
        {
            const auto& curv = opt.getCurrentVertexCurvature(v);
            EXPECT_NEAR(curv.m_kmax, -1.0, 0.1);
            EXPECT_NEAR(curv.m_kmin, -1.0, 0.1);
        }
    }
}

TEST_F(CurvatureOptimizerCylinderTest, Plane)
{
    CurvatureOptimizer opt(mesh);
    opt.initMesh();
    for (pmp::Vertex v : mesh->vertices())
    {
        opt.addTargetVertexCurvature(v, 0.0, 0.0);
    }
    opt.compute();

    for (pmp::Vertex v : mesh->vertices())
    {
        const auto& curv = opt.getCurrentVertexCurvature(v);
        EXPECT_NEAR(curv.m_kmax, 0.0, 0.1);
        EXPECT_NEAR(curv.m_kmin, 0.0, 0.1);
    }
}

TEST_F(CurvatureOptimizerCylinderTest, Sphere)
{
    CurvatureOptimizer opt(mesh);
    opt.initMesh();
    for (pmp::Vertex v : mesh->vertices())
    {
        opt.addTargetVertexCurvature(v, 1.0, 1.0);
    }
    opt.compute();

    for (pmp::Vertex v : mesh->vertices())
    {
        if (!mesh->is_boundary(v))
        {
            const auto& curv = opt.getCurrentVertexCurvature(v);
            EXPECT_NEAR(curv.m_kmax, 1.0, 0.25);
            EXPECT_NEAR(curv.m_kmin, 1.0, 0.25);
        }
    }
}

TEST_F(CurvatureOptimizerCylinderTest, InvSphere)
{
    CurvatureOptimizer opt(mesh);
    opt.initMesh();
    for (pmp::Vertex v : mesh->vertices())
    {
        opt.addTargetVertexCurvature(v, -1.0, -1.0);
    }
    opt.compute();

    for (pmp::Vertex v : mesh->vertices())
    {
        if (!mesh->is_boundary(v))
        {
            const auto& curv = opt.getCurrentVertexCurvature(v);
            EXPECT_NEAR(curv.m_kmax, -1.0, 0.2);
            EXPECT_NEAR(curv.m_kmin, -1.0, 0.2);
        }
    }
}