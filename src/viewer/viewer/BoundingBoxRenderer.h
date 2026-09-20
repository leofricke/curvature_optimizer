#pragma once

#include "CameraState.h"
#include "LineRenderer.h"

#include <pmp/bounding_box.h>

using namespace pmp;

class BoundingBoxRenderer
{
public:
    BoundingBoxRenderer();

    void update();
    void draw(const CameraState& cameraState);

    void addBoundingBox(const BoundingBox& box);
    void clear();

private:
    bool m_dirty;

    std::vector<BoundingBox> m_boxes;
    LineRenderer m_boxRenderer;
};
