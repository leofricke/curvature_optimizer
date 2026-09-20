list (APPEND HEADER_FILES
    implot.h
    implot_internal.h
)

list (APPEND SOURCE_FILES
    implot.cpp
    implot_items.cpp
)

add_library(implot STATIC ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(implot PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(implot PUBLIC imgui)