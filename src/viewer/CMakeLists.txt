list (APPEND HEADER_FILES
    viewer/AdvancedMeshRenderer.h
    viewer/BaseViewer.h
    viewer/BoundingBoxRenderer.h
    viewer/CameraState.h
    viewer/LineShaderSource.h
    viewer/LineRenderer.h
    viewer/PointShaderSource.h
    viewer/PointRenderer.h

    viewer/buffers/IndexBuffer.h
    viewer/buffers/VertexArray.h
    viewer/buffers/VertexBuffer.h

    viewer/color/ColorMap.h
    viewer/color/ColorUtil.h
)

list (APPEND SOURCE_FILES
    viewer/AdvancedMeshRenderer.cpp
    viewer/BaseViewer.cpp
    viewer/BoundingBoxRenderer.cpp
    viewer/LineRenderer.cpp
    viewer/PointRenderer.cpp

    viewer/buffers/IndexBuffer.cpp
    viewer/buffers/VertexArray.cpp
    viewer/buffers/VertexBuffer.cpp

    viewer/color/ColorMap.cpp
    viewer/color/ColorUtil.cpp
)

add_library(viewer STATIC ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(viewer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(viewer PUBLIC pmp_viewers implot)
target_link_libraries(viewer PRIVATE utils io)