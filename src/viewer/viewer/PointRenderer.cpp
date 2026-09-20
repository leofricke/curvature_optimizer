#include "PointRenderer.h"

#include "PointShaderSource.h"

PointRenderer::PointRenderer()
    : m_dirty(true), m_defaultPointColor(0.0f), m_pointSize(8.0f)
{
    m_numInstances = 0;

    m_vertexArray.addBuffer(3, 0, GL_FLOAT);
    m_vertexArray.addBuffer(3, 1, GL_FLOAT);
    m_vertexArray.addBuffer(3, 1, GL_FLOAT);
}

void PointRenderer::update()
{
    m_dirty = true;
}

void PointRenderer::draw(const CameraState& cameraState)
{
    if (!m_shader.is_valid())
    {
        try
        {
            m_shader.source(point_vshader, point_fshader);
        }
        catch (const GLException& e)
        {
            std::cerr << e.what() << std::endl;
            exit(1);
        }
    }

    if (m_dirty)
    {
        m_dirty = false;

        m_numInstances = m_pointPositions.size();
        createGeometry();

        m_vertexArray.bind();

        m_vertexArray.getBuffer(0).upload(m_geometry.data(), m_geometry.size() * 3 * sizeof(float));
        m_vertexArray.getBuffer(1).upload(m_pointPositions.data(), m_numInstances * 3 * sizeof(float));
        m_vertexArray.getBuffer(2).upload(m_pointColors.data(), m_numInstances * 3 * sizeof(float));

        m_vertexArray.unbind();
    }

    if (m_numInstances > 0)
    {
        m_shader.use();
        m_shader.set_uniform("u_mvpMatrix", cameraState.m_modelViewProjectionMatrix);
        m_shader.set_uniform("u_screenWidth", static_cast<float>(cameraState.m_screenWidth));
        m_shader.set_uniform("u_screenHeight", static_cast<float>(cameraState.m_screenHeight));
        m_shader.set_uniform("u_pointSize", m_pointSize);

        glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE);

        m_vertexArray.bind();
        glDrawArraysInstanced(GL_TRIANGLES, 0, m_geometry.size(), m_numInstances);
        m_vertexArray.unbind();

        glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);

        m_shader.disable();
    }
}

const vec3& PointRenderer::getDefaultPointColor() const
{
    return m_defaultPointColor;
}

void PointRenderer::setDefaultPointColor(const vec3& color)
{
    m_defaultPointColor = color;
}

float PointRenderer::getPointSize() const
{
    return m_pointSize;
}

void PointRenderer::setPointSize(float size)
{
    m_pointSize = size;
}

void PointRenderer::clear()
{
    m_pointPositions.clear();
    m_pointColors.clear();
}

void PointRenderer::addPoint(const dvec3& p, const std::optional<Color>& color)
{
    addPoint(static_cast<vec3>(p), color);
}

void PointRenderer::addPoint(const vec3& p, const std::optional<Color>& color)
{
    m_pointPositions.emplace_back(p);
    m_pointColors.emplace_back(color ? *color : m_defaultPointColor);
}

void PointRenderer::addPoints(const std::vector<dvec3>& points)
{
    m_pointPositions.reserve(m_pointPositions.size() + points.size());
    m_pointColors.reserve(m_pointColors.size() + points.size());
    for (const dvec3& p : points)
    {
        addPoint(p);
    }
}

void PointRenderer::createGeometry()
{
    m_geometry.clear();

    const size_t resolution = m_pointSize;
    for (size_t i = 0; i < resolution; ++i)
    {
        const float theta0 = i * 2.0 * std::numbers::pi / resolution;
        const float theta1 = (i + 1) * 2.0 * std::numbers::pi / resolution;

        m_geometry.emplace_back(0.0, 0.0, 0.0);
        m_geometry.emplace_back(0.5 * std::cos(theta0), 0.5 * std::sin(theta0), 0.0);
        m_geometry.emplace_back(0.5 * std::cos(theta1), 0.5 * std::sin(theta1), 0.0);
    }
}