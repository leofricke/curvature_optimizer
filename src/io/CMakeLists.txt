list (APPEND HEADER_FILES
    io/JsonPMP.h
    io/MeshIO.h
    io/MeshSequenceExporter.h
    io/ReadPLY.h
    io/WritePLY.h
)

list (APPEND SOURCE_FILES
    io/ReadPLY.cpp
    io/MeshIO.cpp
    io/MeshSequenceExporter.cpp
    io/WritePLY.cpp
)

add_library(io ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(io PUBLIC utils nlohmann_json happly)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(io PUBLIC OpenMP::OpenMP_CXX)
endif()

if(WIN32)
    set_property(TARGET io PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()