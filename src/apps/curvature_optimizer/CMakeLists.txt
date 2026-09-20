list (APPEND HEADER_FILES
    CurvatureOptimizerViewer.h
)

list (APPEND SOURCE_FILES
    Main.cpp
    CurvatureOptimizerViewer.cpp
)

if(NOT BUILD_HEADLESS)
    add_app(curvature_optimizer_app "${SOURCE_FILES}")
    target_include_directories(curvature_optimizer_app PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(curvature_optimizer_app PRIVATE curvature_optimizer io viewer args)
endif()