#include "LineRenderer.h"

#include "LineShaderSource.h"

LineRenderer::LineRenderer()
    : m_dirty(true), m_lineWidth(2.0f), m_lineAlpha(1.0f), m_defaultLineColor(0.0f), m_useRoundCaps(true)
{
    m_numInstances = 0;

    m_vertexArray.addBuffer(3, 0, GL_FLOAT);
    m_vertexArray.addBuffer(3, 1, GL_FLOAT);
    m_vertexArray.addBuffer(3, 1, GL_FLOAT);
    m_vertexArray.addBuffer(3, 1, GL_FLOAT);
}

void LineRenderer::update()
{
    m_dirty = true;
}

void LineRenderer::draw(const CameraState& cameraState)
{
    if (!m_shader.is_valid())
    {
        try
        {
            m_shader.source(line_vshader, line_fshader);
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

        m_numInstances = m_segmentStart.size();
        createGeometry();

        m_vertexArray.bind();

        m_vertexArray.getBuffer(0).upload(m_geometry.data(), m_geometry.size() * 3 * sizeof(float));
        m_vertexArray.getBuffer(1).upload(m_segmentStart.data(), m_numInstances * 3 * sizeof(float));
        m_vertexArray.getBuffer(2).upload(m_segmentEnd.data(), m_numInstances * 3 * sizeof(float));
        m_vertexArray.getBuffer(3).upload(m_segmentColor.data(), m_numInstances * 3 * sizeof(float));

        m_vertexArray.unbind();
    }

    if (m_numInstances > 0)
    {
        m_shader.use();
        m_shader.set_uniform("u_mvpMatrix", cameraState.m_modelViewProjectionMatrix);
        m_shader.set_uniform("u_screenWidth", static_cast<float>(cameraState.m_screenWidth));
        m_shader.set_uniform("u_screenHeight", static_cast<float>(cameraState.m_screenHeight));
        m_shader.set_uniform("u_lineWidth", m_lineWidth);
        m_shader.set_uniform("u_lineAlpha", m_lineAlpha);

        glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE);

        m_vertexArray.bind();
        glDrawArraysInstanced(GL_TRIANGLES, 0, m_geometry.size(), m_numInstances);
        m_vertexArray.unbind();

        glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);

        m_shader.disable();
    }
}

float LineRenderer::getLineWidth() const
{
    return m_lineWidth;
}

void LineRenderer::setLineWidth(float width)
{
    m_lineWidth = width;
}

float LineRenderer::getLineAlpha() const
{
    return m_lineAlpha;
}

void LineRenderer::setLineAlpha(float alpha)
{
    m_lineAlpha = alpha;
}

const vec3& LineRenderer::getDefaultLineColor() const
{
    return m_defaultLineColor;
}

void LineRenderer::setDefaultLineColor(const vec3& color)
{
    m_defaultLineColor = color;
}

bool LineRenderer::getUseRoundCaps() const
{
    return m_useRoundCaps;
}

void LineRenderer::setUseRoundCaps(bool flag)
{
    m_useRoundCaps = flag;
}

void LineRenderer::addSegment(const dvec3& start, const dvec3& end, const std::optional<Color>& color)
{
    addSegment(static_cast<vec3>(start), static_cast<vec3>(end), color);
}

void LineRenderer::addSegment(const vec3& start, const vec3& end, const std::optional<Color>& color)
{
    m_segmentStart.emplace_back(start);
    m_segmentEnd.emplace_back(end);
    m_segmentColor.emplace_back(color ? *color : m_defaultLineColor);
}

void LineRenderer::clear()
{
    m_segmentStart.clear();
    m_segmentEnd.clear();
    m_segmentColor.clear();
}

void LineRenderer::createGeometry()
{
    m_geometry.clear();

    m_geometry.emplace_back(0.0, -0.5, 0.0);
    m_geometry.emplace_back(0.0, -0.5, 1.0);
    m_geometry.emplace_back(0.0, 0.5, 1.0);
    m_geometry.emplace_back(0.0, -0.5, 0.0);
    m_geometry.emplace_back(0.0, 0.5, 1.0);
    m_geometry.emplace_back(0.0, 0.5, 0.0);

    if (m_useRoundCaps && m_lineWidth > 2.0f)
    {
        const size_t capResolution = static_cast<size_t>(m_lineWidth / 2.0);
        for (size_t i = 0; i < capResolution; ++i)
        {
            const float theta0 = std::numbers::pi / 2.0 + (i * std::numbers::pi) / capResolution;
            const float theta1 = std::numbers::pi / 2.0 + ((i + 1) * std::numbers::pi) / capResolution;

            m_geometry.emplace_back(0.0, 0.0, 0.0);
            m_geometry.emplace_back(0.5 * std::cos(theta0), 0.5 * std::sin(theta0), 0.0);
            m_geometry.emplace_back(0.5 * std::cos(theta1), 0.5 * std::sin(theta1), 0.0);
        }
        for (size_t i = 0; i < capResolution; ++i)
        {
            const float theta0 = 3.0 * std::numbers::pi / 2.0 + (i * std::numbers::pi) / capResolution;
            const float theta1 = 3.0 * std::numbers::pi / 2.0 + ((i + 1) * std::numbers::pi) / capResolution;

            m_geometry.emplace_back(0.0, 0.0, 1.0);
            m_geometry.emplace_back(0.5 * std::cos(theta0), 0.5 * std::sin(theta0), 1.0);
            m_geometry.emplace_back(0.5 * std::cos(theta1), 0.5 * std::sin(theta1), 1.0);
        }
    }
}