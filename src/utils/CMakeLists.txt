list (APPEND HEADER_FILES
    
    utils/Barycentric.h
    utils/EigenMathUtil.h
    utils/InlineMatrix.h
    utils/InlineVector.h
    utils/KoppEigendecomposition.h
    utils/MathUtil.h
    utils/Ray.h
    utils/RobustSqrt.h
    utils/Timer.h
)

list (APPEND SOURCE_FILES
    utils/KoppEigendecomposition.cpp
    utils/MathUtil.cpp
)

add_library(utils ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(utils PUBLIC pmp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(utils PUBLIC OpenMP::OpenMP_CXX)
endif()

if(WIN32)
    set_property(TARGET utils PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()