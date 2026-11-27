list (APPEND HEADER_FILES
    opt/EnergyFunction.h
    opt/GaussNewtonEnergyFunction.h
    opt/GaussNewtonOptimizer.h
    opt/GaussNewtonLambdaModule.h
    opt/LinearSolver.h
    opt/LineSearch.h
    opt/MeshDofs.h
    opt/OptimizerBase.h

    opt/utils/SpMatUtil.h
)

list (APPEND SOURCE_FILES
    opt/GaussNewtonEnergyFunction.cpp
    opt/GaussNewtonOptimizer.cpp
    opt/LinearSolver.cpp
    opt/LineSearch.cpp
    opt/OptimizerBase.cpp
)

add_library(opt ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(opt PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(opt PUBLIC utils Eigen3::Eigen spmat::spmat)

include(FindSolverLibraries)
link_solver_libraries(opt)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(opt PUBLIC OpenMP::OpenMP_CXX)
endif()

if(WIN32)
    set_property(TARGET opt PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()
