#pragma once

#include "CameraState.h"
#include "buffers/VertexArray.h"
#include <pmp/viewers/shader.h>
#include <optional>
#include <vector>

using namespace pmp;

class PointRenderer
{
public:
    PointRenderer();

    void update();
    void draw(const CameraState& cameraState);

    const vec3& getDefaultPointColor() const;
    void setDefaultPointColor(const vec3& color);
    float getPointSize() const;
    void setPointSize(float size);

    void clear();
    void addPoint(const dvec3& p, const std::optional<Color>& color = std::nullopt);
    void addPoint(const vec3& p, const std::optional<Color>& color = std::nullopt);
    void addPoints(const std::vector<dvec3>& points);

private:
    void createGeometry();

private:
    std::vector<vec3> m_geometry;
    std::vector<vec3> m_pointPositions;
    std::vector<vec3> m_pointColors;
    bool m_dirty;

    vec3 m_defaultPointColor;
    float m_pointSize;

    Shader m_shader;
    size_t m_numInstances;
    VertexArray m_vertexArray;
};
