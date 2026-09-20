#include "AdvancedMeshRenderer.h"

#include <pmp/algorithms/utilities.h>
#include "color/ColorMap.h"

AdvancedMeshRenderer::AdvancedMeshRenderer(SurfaceMesh& mesh)
    : m_mesh(mesh), m_dirty(true), m_drawMode(DrawMode::Plain), m_meshColor(0.0), m_meshAlpha(1.0), m_showNormals(false), m_showBoundingBox(false), m_meshRenderer(m_mesh)
{
}

void AdvancedMeshRenderer::update()
{
    m_dirty = true;
}

void AdvancedMeshRenderer::draw(const CameraState& cameraState)
{
    if (m_dirty)
    {
        m_dirty = false;

        if (m_drawMode == DrawMode::VertexColors)
        {
            clearFaceColors();
        }
        else if (m_drawMode == DrawMode::FaceColors)
        {
            clearVertexColors();
        }
        m_meshRenderer.set_use_colors(m_drawMode == DrawMode::VertexColors || m_drawMode == DrawMode::FaceColors);
        m_meshRenderer.update_opengl_buffers();

        const auto normals = m_mesh.get_vertex_property<vec3>("v:normal");
        if (normals)
        {
            const float length = 0.33f * mean_edge_length(m_mesh);

            m_normalRenderer.clear();
            for (Vertex v : m_mesh.vertices())
            {
                const vec3& n = normals[v];
                const vec3& p = m_mesh.position(v);

                m_normalRenderer.addSegment(p, p + length * n);
            }
            m_normalRenderer.update();
        }

        m_boundingBoxRenderer.clear();
        m_boundingBoxRenderer.addBoundingBox(bounds(m_mesh));
        m_boundingBoxRenderer.update();
    }

    m_meshRenderer.set_front_color(m_meshColor);
    m_meshRenderer.set_back_color(m_meshColor);
    m_meshRenderer.set_specular(0.0);
    m_meshRenderer.set_alpha(m_meshAlpha);

    glDepthRange(0.01, 1.0);
    m_meshRenderer.draw(cameraState.m_projectionMatrix, cameraState.m_modelViewMatrix, getPmpDrawMode());
    glDepthRange(0.0, 1.0);

    m_vectorFieldRenderer.draw(cameraState);

    if (m_showNormals)
    {
        m_normalRenderer.draw(cameraState);
    }

    if (m_showBoundingBox)
    {
        m_boundingBoxRenderer.draw(cameraState);
    }
}

void AdvancedMeshRenderer::setupVertexColors(const std::function<Color(Vertex)>& func)
{
    if (m_mesh.is_empty())
        return;
    VertexProperty<Color> colors = m_mesh.vertex_property<Color>("v:color");
    for (Vertex v : m_mesh.vertices())
    {
        colors[v] = func(v);
    }

    m_dirty = true;
}

void AdvancedMeshRenderer::setupVertexColorsScalar(const std::function<double(Vertex)>& func)
{
    if (m_mesh.is_empty())
        return;
    std::vector<double> values;
    values.reserve(m_mesh.n_vertices());
    for (Vertex v : m_mesh.vertices())
    {
        values.emplace_back(func(v));
    }
    std::ranges::sort(values);
    const double vmin = values.front();
    const double vmax = values.back();

    ViridisColorMap colorMap;
    setupVertexColors([&](Vertex v) { return colorMap.value((func(v) - vmin) / (vmax - vmin)); });
}

void AdvancedMeshRenderer::clearVertexColors()
{
    if (!m_mesh.has_vertex_property("v:color"))
        return;

    auto colors = m_mesh.get_vertex_property<Color>("v:color");
    m_mesh.remove_vertex_property(colors);

    m_dirty = true;
}

void AdvancedMeshRenderer::setupFaceColors(const std::function<Color(Face)>& func)
{
    if (m_mesh.is_empty())
        return;
    FaceProperty<Color> colors = m_mesh.face_property<Color>("f:color");
    for (Face f : m_mesh.faces())
    {
        colors[f] = func(f);
    }

    m_dirty = true;
}

void AdvancedMeshRenderer::setupFaceColorsScalar(const std::function<double(Face)>& func)
{
    if (m_mesh.is_empty())
        return;
    std::vector<double> values;
    values.reserve(m_mesh.n_faces());
    for (Face f : m_mesh.faces())
    {
        values.emplace_back(func(f));
    }
    std::ranges::sort(values);
    const double vmin = values.front();
    const double vmax = values.back();

    ViridisColorMap colorMap;
    setupFaceColors([&](Face f) { return colorMap.value((func(f) - vmin) / (vmax - vmin)); });
}

void AdvancedMeshRenderer::clearFaceColors()
{
    if (!m_mesh.has_face_property("f:color"))
        return;

    FaceProperty<Color> colors = m_mesh.get_face_property<Color>("f:color");
    m_mesh.remove_face_property(colors);

    m_dirty = true;
}

void AdvancedMeshRenderer::setupVectorField(const std::function<vec3(Vertex)>& func, Color c, bool sym)
{
    if (m_mesh.is_empty())
        return;
    m_vectorFieldRenderer.clear();
    const float refLength = 0.33f * mean_edge_length(m_mesh);

    std::vector<std::pair<vec3, vec3>> field;
    field.reserve(m_mesh.n_vertices());
    float maxLength = 0.0;
    for (Vertex v : m_mesh.vertices())
    {
        const vec3 v1 = func(v);
        field.emplace_back(m_mesh.position(v), v1);

        maxLength = std::max(maxLength, norm(v1));
    }

    for (const auto& [point, vec] : field)
    {
        const Scalar l = refLength * (norm(vec) / maxLength);
        const vec3 d = l * normalize(vec);
        if (sym)
        {
            m_vectorFieldRenderer.addSegment(point - d, point + d, c);
        }
        else
        {
            m_vectorFieldRenderer.addSegment(point, point + d, c);
        }
    }
    m_vectorFieldRenderer.update();
}

void AdvancedMeshRenderer::setupVectorCrossField(const std::function<std::pair<vec3, vec3>(Vertex)>& func, Color c1, Color c2)
{
    if (m_mesh.is_empty())
        return;
    m_vectorFieldRenderer.clear();
    const float refLength = 0.33f * mean_edge_length(m_mesh);

    std::vector<std::tuple<vec3, vec3, vec3>> field;
    field.reserve(m_mesh.n_vertices());
    float maxLength = 0.0;
    for (Vertex v : m_mesh.vertices())
    {
        const auto [v1, v2] = func(v);
        field.emplace_back(m_mesh.position(v), v1, v2);

        maxLength = std::max({ maxLength, norm(v1), norm(v2) });
    }

    for (const auto& [point, v1, v2] : field)
    {
        const Scalar l1 = refLength * (norm(v1) / maxLength);
        const Scalar l2 = refLength * (norm(v2) / maxLength);
        const vec3 d1 = l1 * normalize(v1);
        const vec3 d2 = l2 * normalize(v2);
        m_vectorFieldRenderer.addSegment(point - d1, point + d1, c1);
        m_vectorFieldRenderer.addSegment(point - d2, point + d2, c2);
    }
    m_vectorFieldRenderer.update();
}

void AdvancedMeshRenderer::clearVectorField()
{
    m_vectorFieldRenderer.clear();
    m_vectorFieldRenderer.update();
}

AdvancedMeshRenderer::DrawMode AdvancedMeshRenderer::getDrawMode() const
{
    return m_drawMode;
}

void AdvancedMeshRenderer::setDrawMode(DrawMode mode)
{
    m_drawMode = mode;
    m_dirty = true;
}

bool AdvancedMeshRenderer::getShowNormals() const
{
    return m_showNormals;
}

void AdvancedMeshRenderer::setShowNormals(bool show)
{
    m_showNormals = show;
}

bool AdvancedMeshRenderer::getShowBoundingBox() const
{
    return m_showBoundingBox;
}

void AdvancedMeshRenderer::setShowBoundingBox(bool show)
{
    m_showBoundingBox = show;
}

const vec3& AdvancedMeshRenderer::getMeshColor() const
{
    return m_meshColor;
}

void AdvancedMeshRenderer::setMeshColor(const vec3& color)
{
    m_meshColor = color;
}

float AdvancedMeshRenderer::getMeshAlpha() const
{
    return m_meshAlpha;
}

void AdvancedMeshRenderer::setMeshAlpha(float alpha)
{
    m_meshAlpha = alpha;
}

float AdvancedMeshRenderer::getCreaseAngle() const
{
    return m_meshRenderer.crease_angle();
}

void AdvancedMeshRenderer::setCreaseAngle(float angle)
{
    m_meshRenderer.set_crease_angle(angle);
}

std::string AdvancedMeshRenderer::getPmpDrawMode() const
{
    switch (m_drawMode)
    {
    case DrawMode::Plain:
        return "Smooth Shading";
        break;
    case DrawMode::Faces:
        return "Hidden Line";
        break;
    case DrawMode::VertexColors:
        return "Smooth Shading";
        break;
    case DrawMode::FaceColors:
        return "Smooth Shading";
        break;
    default:
        throw std::runtime_error("Unknown draw mode");
    }
}
