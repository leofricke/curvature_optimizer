#include "BaseViewer.h"

#include "io/JsonPMP.h"

#include "implot.h"
#include <pmp/algorithms/utilities.h>

BaseViewer::BaseViewer(const char* title, int width, int height, bool showgui)
    : TrackballViewer(title, width, height, showgui)
{
    ImPlot::CreateContext();

    m_mesh = std::make_shared<SurfaceMesh>();
    m_meshRenderer = std::make_unique<AdvancedMeshRenderer>(*m_mesh);

    clear_draw_modes();
}

void BaseViewer::process_imgui()
{
    drawLabels();
}

void BaseViewer::drawLabels() const
{
    const ImColor& textColor = ImGui::GetStyle().Colors[ImGuiCol_Text];
    const ImColor& bgColor = ImGui::GetStyle().Colors[ImGuiCol_ButtonActive];
    ImDrawList* drawList = ImGui::GetBackgroundDrawList();

    for (size_t i = 0; i < m_labels.size(); ++i)
    {
        const auto& label = m_labels[i];
        double x, y, depth;
        worldToPixelCoords(label.m_position, x, y, depth);

        if (depth >= 0.0)
        {
            const ImVec2 textSize = ImGui::CalcTextSize(label.m_text.c_str());
            drawList->AddRectFilled(ImVec2(x, y), ImVec2(x + textSize[0] + 5, y + textSize[1] + 4), bgColor, 4);
            drawList->AddText(ImVec2(x + 3, y + 2), textColor, label.m_text.c_str());
        }
    }
}

void BaseViewer::addLabel(const vec3& position, const std::string& text)
{
    m_labels.emplace_back(position, text);
}

void BaseViewer::clearLabels()
{
    m_labels.clear();
}

void BaseViewer::worldToPixelCoords(const vec3& point, double& x, double& y, double& depth) const
{
    vec4 p(point[0], point[1], point[2], 1.0f);
    p = projection_matrix_ * modelview_matrix_ * p;
    depth = p[2];
    p /= p[3];

    x = 0.5 * p[0] * width() + 0.5 * width();
    y = -0.5 * p[1] * height() + 0.5 * height();
}

Ray BaseViewer::createRay(float x, float y) const
{
    const vec2 pos(2.0f * x / width() - 1.0f, 1.0f - 2.0f * y / height());
    const mat4 inv = inverse(projection_matrix_ * modelview_matrix_);

    vec4 nearPoint = inv * vec4(pos[0], pos[1], -1.0f, 1.0f);
    vec4 farPoint = inv * vec4(pos[0], pos[1], 1.0f, 1.0f);
    nearPoint /= nearPoint[3];
    farPoint /= farPoint[3];

    const vec3 origin = vec3(nearPoint[0], nearPoint[1], nearPoint[2]);
    const vec3 dir = vec3(farPoint[0] - nearPoint[0], farPoint[1] - nearPoint[1], farPoint[2] - nearPoint[2]);

    return Ray(origin, normalize(dir));
}

std::pair<Face, vec3> BaseViewer::pickFace(float x, float y, const std::function<bool(Face)>& accept) const
{
    const Ray ray = createRay(x, y);

    Face faceHit;
    Scalar tHit(std::numeric_limits<Scalar>::max());
    vec3 posHit(std::numeric_limits<Scalar>::max());
    for (Face f : m_mesh->faces())
    {
        if (m_mesh->valence(f) != 3)
            continue;
        if (accept && !accept(f))
            continue;

        Halfedge h = m_mesh->halfedge(f);
        const Point a = m_mesh->position(m_mesh->to_vertex(h));
        h = m_mesh->next_halfedge(h);
        const Point b = m_mesh->position(m_mesh->to_vertex(h));
        h = m_mesh->next_halfedge(h);
        const Point c = m_mesh->position(m_mesh->to_vertex(h));

        Scalar t;
        vec3 pointHit;
        const bool hit = ray.intersectTriangle(a, b, c, t, pointHit);
        if (hit && t < tHit)
        {
            tHit = t;
            posHit = ray.evaluate(t);
            faceHit = f;
        }
    }

    return { faceHit, posHit };
}

BarycentricPoint BaseViewer::pickFaceBarycentric(float x, float y, const std::function<bool(Face)>& accept) const
{
    const auto [face, pos] = pickFace(x, y, accept);
    if (face.is_valid())
    {
        return BarycentricPoint::compute(*m_mesh, m_mesh->halfedge(face), pos);
    }
    return BarycentricPoint();
}

Vertex BaseViewer::pickVertex(float x, float y, const std::function<bool(Vertex)>& accept, float pickingRadius) const
{
    const Scalar mean = mean_edge_length(*m_mesh);
    if (pickingRadius <= 0.0f) // automatic
    {
        pickingRadius = mean * 0.3f;
    }

    const Ray ray = createRay(x, y);

    const auto [face, pos] = pickFace(x, y);
    const Scalar tMax = distance(ray.m_origin, pos);

    Vertex vertexHit;
    Scalar tHit(std::numeric_limits<Scalar>::max());
    for (Vertex v : m_mesh->vertices())
    {
        if (accept && !accept(v))
            continue;

        Scalar t;
        vec3 pointHit;
        const bool hit = ray.intersectSphere(m_mesh->position(v), pickingRadius, t, pointHit);

        if (hit && t < tHit && t <= tMax)
        {
            vertexHit = v;
            tHit = t;
        }
    }
    return vertexHit;
}

std::string BaseViewer::getCameraStateAsJson() const
{
    const json j = {
        { "center", center_ },
        { "radius", radius_ },
        { "modelViewMatrix", modelview_matrix_ },
    };

    return j.dump();
}

void BaseViewer::setCameraStateAsJson(const std::string& jsonString)
{
    try
    {
        const json j = json::parse(jsonString);

        center_ = j.at("center");
        radius_ = j.at("radius");
        modelview_matrix_ = j.at("modelViewMatrix");
    }
    catch (const std::exception&)
    {
        std::cout << "ERROR: Failed to parse camera state" << std::endl;
    }
}
