#pragma once

#include "CameraState.h"
#include "LineRenderer.h"
#include "BoundingBoxRenderer.h"

#include <pmp/surface_mesh.h>
#include <pmp/viewers/renderer.h>

using namespace pmp;

class AdvancedMeshRenderer
{
public:
    enum class DrawMode
    {
        Plain = 0,
        Faces = 1,
        VertexColors = 2,
        FaceColors = 3,
    };

public:
    explicit AdvancedMeshRenderer(SurfaceMesh& mesh);

    void update();
    void draw(const CameraState& cameraState);

    void setupVertexColors(const std::function<Color(Vertex)>& func);
    void setupVertexColorsScalar(const std::function<double(Vertex)>& func);
    void clearVertexColors();

    void setupFaceColors(const std::function<Color(Face)>& func);
    void setupFaceColorsScalar(const std::function<double(Face)>& func);
    void clearFaceColors();

    void setupVectorField(const std::function<vec3(Vertex)>& func, Color c = Color(0, 0, 0), bool sym = false);
    void setupVectorCrossField(const std::function<std::pair<vec3, vec3>(Vertex)>& func, Color c1 = Color(1, 0, 0), Color c2 = Color(0, 0, 1));
    void clearVectorField();

    DrawMode getDrawMode() const;
    void setDrawMode(DrawMode mode);

    bool getShowNormals() const;
    void setShowNormals(bool show);

    bool getShowBoundingBox() const;
    void setShowBoundingBox(bool show);

    const vec3& getMeshColor() const;
    void setMeshColor(const vec3& color);
    float getMeshAlpha() const;
    void setMeshAlpha(float alpha);

    float getCreaseAngle() const;
    void setCreaseAngle(float angle);

private:
    std::string getPmpDrawMode() const;

private:
    SurfaceMesh& m_mesh;
    bool m_dirty;

    DrawMode m_drawMode;
    vec3 m_meshColor;
    float m_meshAlpha;
    bool m_showNormals;
    bool m_showBoundingBox;

    Renderer m_meshRenderer;

    LineRenderer m_vectorFieldRenderer;
    LineRenderer m_normalRenderer;
    BoundingBoxRenderer m_boundingBoxRenderer;
};
