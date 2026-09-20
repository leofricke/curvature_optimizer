#pragma once

#include "viewer/BaseViewer.h"

#include "curvature_optimizer/CurvatureOptimizer.h"
#include "viewer/AdvancedMeshRenderer.h"
#include "io/MeshSequenceExporter.h"

#include <atomic>
#include <filesystem>

enum class DrawMode
{
    Plain = 0,
    Faces = 1,
    Curvature = 2
};
static std::array<std::string, 3> s_drawModeNames = { "Plain", "Faces", "Curvature" };

enum class CurvatureDrawMode
{
    Mean = 0,
    Gaussian = 1,
    MaxAbs = 2
};
static std::array<std::string, 3> s_curvatureDrawModeNames = { "Mean", "Gaussian", "MaxAbs" };

enum class CurvatureEditMode
{
    Global = 0,
    Scale = 1,
    FeatureEnhance = 2,
    Clamp = 3
};
static std::array<std::string, 4> s_curvatureEditModeNames = { "Global", "Scale", "FeatureEnhance", "Clamp" };

enum class BoundaryCurvatureEditMode
{
    Global = 0,
    Scale = 1
};
static std::array<std::string, 2> s_boundaryCurvatureEditModeNames = { "Global", "Scale" };

class CurvatureOptimizerViewer : public BaseViewer
{
public:
    CurvatureOptimizerViewer(const char* title, int width, int height, bool hideGuiOnStart = false, bool collapseGuiOnStart = false);

    void loadMesh(std::filesystem::path meshFile);
    void updateMesh(bool reset = true);
    void resetDefaultState();
    void updateMeshDrawMode();

    void saveCurrentCurvatures();
    void startOptimizer(bool useSavedCurvs);

    void drop(int count, const char** paths) override;
    void draw(const std::string&) override;
    void process_imgui() override;
    void mouse(int button, int action, int mods) override;
    void keyboard(int key, int code, int action, int mods) override;

protected:
    std::filesystem::path m_meshFile;

    CurvatureOptimizer m_optimizer;
    std::atomic_bool m_running;
    std::atomic_bool m_meshDirty;
    std::map<std::string, std::vector<double>> m_plots;
    MeshSequenceExporter m_sequenceExporter;
    bool m_writeSequence;

    CurvatureEditMode m_curvatureEditMode = CurvatureEditMode::Global;
    BoundaryCurvatureEditMode m_boundaryCurvatureEditMode = BoundaryCurvatureEditMode::Global;
    bool m_lockBoundary = false;
    bool m_lockTwoRing = false;
    bool m_optimizeBoundary = true;
    double m_pNormFairing = 2.0;
    bool m_useBoundaryCurvs = false;
    double m_globalMaxCurv = 0.0;
    double m_globalMinCurv = 0.0;
    double m_globalBoundaryCurv = 0.0;
    double m_scaleMaxCurv = 1.0;
    double m_scaleMinCurv = 1.0;
    double m_scaleBoundaryCurv = 1.0;
    double m_featureEnhanceScale = 1.0;
    double m_clampMin = 0.0;
    double m_clampMax = 0.0;

    bool m_collapseGuiOnStart;
    DrawMode m_drawMode;
    CurvatureDrawMode m_curvatureDrawMode;
    bool m_showCurvatureDirections;
};