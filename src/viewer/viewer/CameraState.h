#pragma once

#include <pmp/surface_mesh.h>

struct CameraState
{
    pmp::mat4 m_projectionMatrix;
    pmp::mat4 m_modelViewMatrix;
    pmp::mat4 m_modelViewProjectionMatrix;
    int m_screenWidth;
    int m_screenHeight;
};