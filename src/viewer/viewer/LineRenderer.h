#pragma once

#include "CameraState.h"
#include "buffers/VertexArray.h"
#include <pmp/viewers/shader.h>
#include <optional>
#include <vector>

using namespace pmp;

class LineRenderer
{
public:
    LineRenderer();

    void update();
    void draw(const CameraState& cameraState);

    float getLineWidth() const;
    void setLineWidth(float width);
    float getLineAlpha() const;
    void setLineAlpha(float alpha);
    const vec3& getDefaultLineColor() const;
    void setDefaultLineColor(const vec3& color);
    bool getUseRoundCaps() const;
    void setUseRoundCaps(bool flag);

    void addSegment(const dvec3& start, const dvec3& end, const std::optional<Color>& color = std::nullopt);
    void addSegment(const vec3& start, const vec3& end, const std::optional<Color>& color = std::nullopt);
    void clear();

private:
    void createGeometry();

private:
    std::vector<vec3> m_geometry;
    std::vector<vec3> m_segmentStart;
    std::vector<vec3> m_segmentEnd;
    std::vector<vec3> m_segmentColor;
    bool m_dirty;
    float m_lineWidth;
    float m_lineAlpha;
    vec3 m_defaultLineColor;
    bool m_useRoundCaps;

    Shader m_shader;
    size_t m_numInstances;
    VertexArray m_vertexArray;
};
