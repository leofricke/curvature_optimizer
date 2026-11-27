list (APPEND HEADER_FILES
    curvature_optimizer/CurvatureTensorVertexCurvature.h
    curvature_optimizer/IVertexCurvatureOperator.h
    curvature_optimizer/CurvatureOptimizer.h
    curvature_optimizer/CurvatureOptimizerWorkspace.h
    curvature_optimizer/MeanGaussianVertexCurvature.h
    curvature_optimizer/PrincipalCurvatures.h
    curvature_optimizer/TargetData.h
    curvature_optimizer/Tolerances.h
    curvature_optimizer/VirtualSubdivision.h
    
    curvature_optimizer/energies/AreaEnergy.h
    curvature_optimizer/energies/BoundaryCurvatureEnergy.h
    curvature_optimizer/energies/ConformalEnergy.h
    curvature_optimizer/energies/EnergyUtil.h
    curvature_optimizer/energies/CurvatureEnergy.h
    curvature_optimizer/energies/CurvatureVariationEnergy.h
    curvature_optimizer/energies/PolygonPlanarityEnergy.h
    curvature_optimizer/energies/PositionEnergy.h
    curvature_optimizer/energies/RegularizationEnergy.h
    curvature_optimizer/energies/WillmoreEnergy.h
    
    curvature_optimizer/quantities/BoundaryCurvatureDerivatives.h
    curvature_optimizer/quantities/CommonDerivatives.h
    curvature_optimizer/quantities/DiffRules.h
    curvature_optimizer/quantities/DihedralAngleDerivatives.h
    curvature_optimizer/quantities/DualAreaDerivatives.h
    curvature_optimizer/quantities/FaceAreaDerivatives.h
    curvature_optimizer/quantities/GaussianCurvatureDerivatives.h
    curvature_optimizer/quantities/InteriorAngleDerivatives.h
    curvature_optimizer/quantities/LocalGradient.h
    curvature_optimizer/quantities/MeanCurvatureDerivatives.h
    curvature_optimizer/quantities/ScalarFieldGradientDerivatives.h
)

list (APPEND SOURCE_FILES
    curvature_optimizer/CurvatureTensorVertexCurvature.cpp
    curvature_optimizer/CurvatureOptimizer.cpp
    curvature_optimizer/CurvatureOptimizerWorkspace.cpp
    curvature_optimizer/MeanGaussianVertexCurvature.cpp
    curvature_optimizer/TargetData.cpp
    curvature_optimizer/VirtualSubdivision.cpp
    
    curvature_optimizer/energies/AreaEnergy.cpp
    curvature_optimizer/energies/BoundaryCurvatureEnergy.cpp
    curvature_optimizer/energies/ConformalEnergy.cpp
    curvature_optimizer/energies/CurvatureEnergy.cpp
    curvature_optimizer/energies/CurvatureVariationEnergy.cpp
    curvature_optimizer/energies/PolygonPlanarityEnergy.cpp
    curvature_optimizer/energies/PositionEnergy.cpp
    curvature_optimizer/energies/RegularizationEnergy.cpp
    curvature_optimizer/energies/WillmoreEnergy.cpp
    
    curvature_optimizer/quantities/BoundaryCurvatureDerivatives.cpp
    curvature_optimizer/quantities/CommonDerivatives.cpp
    curvature_optimizer/quantities/DihedralAngleDerivatives.cpp
    curvature_optimizer/quantities/DualAreaDerivatives.cpp
    curvature_optimizer/quantities/FaceAreaDerivatives.cpp
    curvature_optimizer/quantities/GaussianCurvatureDerivatives.cpp
    curvature_optimizer/quantities/InteriorAngleDerivatives.cpp
    curvature_optimizer/quantities/MeanCurvatureDerivatives.cpp
    curvature_optimizer/quantities/ScalarFieldGradientDerivatives.cpp
)

add_library(curvature_optimizer ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(curvature_optimizer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(curvature_optimizer PUBLIC utils opt)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(curvature_optimizer PUBLIC OpenMP::OpenMP_CXX)
endif()

if(WIN32)
    set_property(TARGET curvature_optimizer PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()