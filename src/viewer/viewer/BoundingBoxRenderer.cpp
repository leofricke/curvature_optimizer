#include "BoundingBoxRenderer.h"

BoundingBoxRenderer::BoundingBoxRenderer()
    : m_dirty(true)
{
}

void BoundingBoxRenderer::update()
{
    m_dirty = true;
}

void BoundingBoxRenderer::draw(const CameraState& cameraState)
{
    if (m_dirty)
    {
        m_dirty = false;

        m_boxRenderer.clear();
        for (const BoundingBox& box : m_boxes)
        {
            const vec3& min = box.min();
            const vec3& max = box.max();

            const vec3 v0 = min;
            const vec3 v1 = vec3(max[0], min[1], min[2]);
            const vec3 v2 = vec3(min[0], max[1], min[2]);
            const vec3 v3 = vec3(max[0], max[1], min[2]);

            const vec3 v4 = vec3(min[0], min[1], max[2]);
            const vec3 v5 = vec3(max[0], min[1], max[2]);
            const vec3 v6 = vec3(min[0], max[1], max[2]);
            const vec3 v7 = max;

            m_boxRenderer.addSegment(v0, v1);
            m_boxRenderer.addSegment(v1, v3);
            m_boxRenderer.addSegment(v3, v2);
            m_boxRenderer.addSegment(v2, v0);

            m_boxRenderer.addSegment(v4, v5);
            m_boxRenderer.addSegment(v5, v7);
            m_boxRenderer.addSegment(v7, v6);
            m_boxRenderer.addSegment(v6, v4);

            m_boxRenderer.addSegment(v0, v4);
            m_boxRenderer.addSegment(v1, v5);
            m_boxRenderer.addSegment(v2, v6);
            m_boxRenderer.addSegment(v3, v7);
        }
        m_boxRenderer.update();
    }

    m_boxRenderer.draw(cameraState);
}

void BoundingBoxRenderer::addBoundingBox(const BoundingBox& box)
{
    m_dirty = true;
    m_boxes.push_back(box);
}

void BoundingBoxRenderer::clear()
{
    m_dirty = true;
    m_boxes.clear();
}
