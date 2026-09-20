#pragma once

#include "utils/Barycentric.h"
#include "utils/Ray.h"
#include "AdvancedMeshRenderer.h"

#include <pmp/surface_mesh.h>
#include <pmp/viewers/trackball_viewer.h>

using namespace pmp;

class BaseViewer : public TrackballViewer
{
private:
    struct Label3D
    {
        vec3 m_position;
        std::string m_text;
    };

public:
    BaseViewer(const char* title, int width, int height, bool showgui = true);

protected:
    virtual void process_imgui() override;

public:
    void drawLabels() const;
    void addLabel(const vec3& position, const std::string& text);
    void clearLabels();

    void worldToPixelCoords(const vec3& point, double& x, double& y, double& depth) const;
    Ray createRay(float x, float y) const;
    std::pair<Face, vec3> pickFace(float x, float y, const std::function<bool(Face)>& accept = nullptr) const;
    BarycentricPoint pickFaceBarycentric(float x, float y, const std::function<bool(Face)>& accept = nullptr) const;
    Vertex pickVertex(float x, float y, const std::function<bool(Vertex)>& accept = nullptr, float pickingRadius = -1) const;

    std::string getCameraStateAsJson() const;
    void setCameraStateAsJson(const std::string& jsonString);

protected:
    std::shared_ptr<SurfaceMesh> m_mesh;
    std::unique_ptr<AdvancedMeshRenderer> m_meshRenderer;

    std::vector<Label3D> m_labels;
};
