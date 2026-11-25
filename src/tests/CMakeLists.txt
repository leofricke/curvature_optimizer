list (APPEND HEADER_FILES
    TestFilesUtil.h
    DerivativeTestUtil.h
    
    opt/PowellFunction.h
    opt/RosenbrockFunction.h
    opt/SimpleQuadraticFunction.h
)

list (APPEND SOURCE_FILES
    curvature_optimizer/TestCurvatureOptimizer.cpp
    curvature_optimizer/TestCurvatureTensorVertexCurvature.cpp
    curvature_optimizer/TestGaussBonnetTheorem.cpp
    curvature_optimizer/TestMeanGaussianVertexCurvature.cpp
    curvature_optimizer/TestVirtualSubdivision.cpp
    
    curvature_optimizer/energies/TestAreaEnergy.cpp
    curvature_optimizer/energies/TestBoundaryCurvatureEnergy.cpp
    curvature_optimizer/energies/TestConformalEnergy.cpp
    curvature_optimizer/energies/TestCurvatureEnergy.cpp
    curvature_optimizer/energies/TestCurvatureVariationEnergy.cpp
    curvature_optimizer/energies/TestPolygonPlanarityEnergy.cpp
    curvature_optimizer/energies/TestPositionEnergy.cpp
    curvature_optimizer/energies/TestRegularizationEnergy.cpp
    
    curvature_optimizer/quantities/TestBoundaryCurvatureDerivatives.cpp
    curvature_optimizer/quantities/TestCommonDerivatives.cpp
    curvature_optimizer/quantities/TestDihedralAngleDerivatives.cpp
    curvature_optimizer/quantities/TestDualAreaDerivatives.cpp
    curvature_optimizer/quantities/TestFaceAreaDerivatives.cpp
    curvature_optimizer/quantities/TestGaussianCurvatureDerivatives.cpp
    curvature_optimizer/quantities/TestInteriorAngleDerivatives.cpp
    curvature_optimizer/quantities/TestMeanCurvatureDerivatives.cpp
    curvature_optimizer/quantities/TestScalarFieldGradientDerivatives.cpp
    
    opt/TestGaussNewtonOptimizer.cpp
    opt/TestGaussNewtonSystem.cpp
    
    utils/TestKoppEigendecomposition.cpp
    utils/TestInlineMatrix.cpp
    utils/TestInlineVector.cpp
)

enable_testing()

add_executable(test_runner ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(test_runner PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(test_runner utils io opt curvature_optimizer)
target_link_libraries(test_runner GTest::gtest_main)
target_compile_definitions(test_runner PRIVATE TEST_FILES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/testfiles")

include(GoogleTest)
gtest_discover_tests(test_runner)