#include "CurvatureOptimizerViewer.h"

#include "io/MeshIO.h"
#include "io/WritePLY.h"

#include <pmp/algorithms/normals.h>
#include <pmp/algorithms/utilities.h>
#include "implot.h"

#include <thread>
#include <format>

const vec3 s_meshColor(0.6, 0.6, 0.6);

CurvatureOptimizerViewer::CurvatureOptimizerViewer(const char* title, int width, int height, bool hideGuiOnStart, bool collapseGuiOnStart)
    : BaseViewer(title, width, height, !hideGuiOnStart),
      m_optimizer(m_mesh),
      m_running(false),
      m_meshDirty(false),
      m_sequenceExporter("mesh_sequence"),
      m_writeSequence(false),
      m_collapseGuiOnStart(collapseGuiOnStart),
      m_drawMode(DrawMode::Plain),
      m_curvatureDrawMode(CurvatureDrawMode::MaxAbs),
      m_showCurvatureDirections(false)
{
    m_optimizer.setVerbose(true);

    m_meshRenderer->setMeshColor(s_meshColor);
    m_meshRenderer->setCreaseAngle(0.0);

    resetDefaultState();

    add_help_item("Backspace", "Reload mesh", 3);
}

void CurvatureOptimizerViewer::loadMesh(std::filesystem::path meshFile)
{
    if (m_running)
    {
        return;
    }

    m_mesh->clear();
    try
    {
        readMesh(meshFile, *m_mesh);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }
    if (m_mesh->is_empty())
    {
        m_meshFile.clear();
    }
    else
    {
        m_meshFile = meshFile;
        std::cout << "Loaded " << meshFile.string() << ": " << m_mesh->n_vertices() << " vertices, " << m_mesh->n_faces() << " faces" << std::endl;
    }

    // update scene center and bounds
    const BoundingBox bb = bounds(*m_mesh);
    set_scene((vec3) bb.center(), 0.5f * bb.size());

    m_plots.clear();

    m_optimizer.initMesh();
    m_optimizer.clearTargets();
    updateMesh();
}

void CurvatureOptimizerViewer::updateMesh(bool reset)
{
    if (reset)
    {
        resetDefaultState();

        // update scene center and radius, but don't update camera view
        const BoundingBox bb = bounds(*m_mesh);
        center_ = (vec3) bb.center();
        radius_ = 2.0f * bb.size();
    }

    vertex_normals(*m_mesh);

    updateMeshDrawMode();

    m_meshRenderer->update();
    clearLabels();
}

void CurvatureOptimizerViewer::resetDefaultState()
{
    m_meshRenderer->setMeshAlpha(1.0);
    m_meshRenderer->setCreaseAngle(0.0);
}

void CurvatureOptimizerViewer::updateMeshDrawMode()
{
    if (m_drawMode == DrawMode::Plain)
    {
        m_meshRenderer->setDrawMode(AdvancedMeshRenderer::DrawMode::Plain);
    }
    else if (m_drawMode == DrawMode::Faces)
    {
        m_meshRenderer->setDrawMode(AdvancedMeshRenderer::DrawMode::Faces);
    }
    else if (m_drawMode == DrawMode::Curvature)
    {
        if (m_curvatureDrawMode == CurvatureDrawMode::Mean)
        {
            m_meshRenderer->setupVertexColorsScalar([&](Vertex v) { return m_optimizer.getCurrentVertexCurvature(v).meanCurvature(); });
        }
        else if (m_curvatureDrawMode == CurvatureDrawMode::Gaussian)
        {
            m_meshRenderer->setupVertexColorsScalar([&](Vertex v) { return m_optimizer.getCurrentVertexCurvature(v).gaussianCurvature(); });
        }
        else if (m_curvatureDrawMode == CurvatureDrawMode::MaxAbs)
        {
            m_meshRenderer->setupVertexColorsScalar([&](Vertex v) { return m_optimizer.getCurrentVertexCurvature(v).maxAbsCurvature(); });
        }
        m_meshRenderer->setDrawMode(AdvancedMeshRenderer::DrawMode::VertexColors);
    }

    if (m_showCurvatureDirections)
    {
        m_meshRenderer->setupVectorCrossField(
            [&](Vertex v)
            {
                const auto& curv = m_optimizer.getCurrentVertexCurvature(v);
                return std::make_pair<vec3, vec3>(curv.m_maxDir, curv.m_minDir);
            });
    }
    else
    {
        m_meshRenderer->clearVectorField();
    }
}

void CurvatureOptimizerViewer::saveCurrentCurvatures()
{
    VertexProperty<Eigen::Vector2d> vertexCurvs = m_mesh->vertex_property<Eigen::Vector2d>("v:savedCurvs", Eigen::Vector2d::Zero());
    VertexProperty<double> boundaryCurvs = m_mesh->vertex_property<double>("v:savedBoundaryCurvs", 0.0);

    for (Vertex v : m_mesh->vertices())
    {
        const auto& curv = m_optimizer.getCurrentVertexCurvature(v);
        vertexCurvs[v] = { curv.m_kmax, curv.m_kmin };
        if (m_mesh->is_boundary(v))
        {
            boundaryCurvs[v] = m_optimizer.getCurrentBoundaryCurvature(v);
        }
    }
}

void CurvatureOptimizerViewer::startOptimizer(bool useSavedCurvs)
{
    if (m_running)
    {
        return;
    }

    m_optimizer.initMesh();
    m_optimizer.clearTargets();

    if (m_writeSequence)
    {
        m_sequenceExporter.reset();
    }
    m_plots.clear();

    // setup callback
    m_optimizer.setMeshUpdateCallback(
        [&](const CurvatureOptimizer::Progress& progress)
        {
            m_plots["Energy"].emplace_back(progress.m_energy);
            for (const auto& [name, energy] : progress.m_blockEnergies)
            {
                m_plots[name].emplace_back(energy);
            }
            if (m_optimizer.getOptimizerMode() == CurvatureOptimizer::OptimizerMode::TargetCurvature)
            {
                m_plots["Curvature Rel Error"].emplace_back(m_optimizer.getCurvatureRelativeError());
            }

            if (m_writeSequence)
            {
                std::vector<double> kmax;
                std::vector<double> kmin;
                for (Vertex v : m_mesh->vertices())
                {
                    const auto& curv = m_optimizer.getCurrentVertexCurvature(v);
                    kmax.push_back(curv.m_kmax);
                    kmin.push_back(curv.m_kmin);
                }
                m_sequenceExporter.appendMesh(*m_mesh, { { "kmax", kmax }, { "kmin", kmin } });
            }

            m_meshDirty = true;
            m_meshDirty.wait(true);
        });

    // setup locked vertices
    std::vector<Vertex> lockedVertices;
    if (m_lockBoundary)
    {
        for (Vertex v : m_mesh->vertices())
        {
            if (m_mesh->is_boundary(v))
            {
                lockedVertices.emplace_back(v);
                if (m_lockTwoRing)
                {
                    for (Vertex vv : m_mesh->vertices(v))
                    {
                        if (!m_mesh->is_boundary(vv) && std::find(lockedVertices.begin(), lockedVertices.end(), vv) == lockedVertices.end())
                        {
                            lockedVertices.emplace_back(vv);
                        }
                    }
                }
            }
        }
    }

    // setup curvature targets
    if (useSavedCurvs)
    {
        if (!m_mesh->has_vertex_property("v:savedCurvs"))
        {
            return;
        }
        VertexProperty<Eigen::Vector2d> vertexCurvs = m_mesh->vertex_property<Eigen::Vector2d>("v:savedCurvs", Eigen::Vector2d::Zero());
        for (Vertex v : m_mesh->vertices())
        {
            if (!m_optimizeBoundary && m_mesh->is_boundary(v))
            {
                continue;
            }
            m_optimizer.addTargetVertexCurvature(v, vertexCurvs[v][0], vertexCurvs[v][1]);
        }
    }
    else
    {
        for (Vertex v : m_mesh->vertices())
        {
            if (!m_optimizeBoundary && m_mesh->is_boundary(v))
            {
                continue;
            }
            const auto& initialCurv = m_optimizer.getCurrentVertexCurvature(v);
            if (m_curvatureEditMode == CurvatureEditMode::Global)
            {
                m_optimizer.addTargetVertexCurvature(v, m_globalMaxCurv, m_globalMinCurv);
            }
            else if (m_curvatureEditMode == CurvatureEditMode::Scale)
            {
                m_optimizer.addTargetVertexCurvature(v, m_scaleMaxCurv * initialCurv.m_kmax, m_scaleMinCurv * initialCurv.m_kmin);
            }
            else if (m_curvatureEditMode == CurvatureEditMode::FeatureEnhance)
            {
                double scaleMax = 1.0;
                double scaleMin = 1.0;
                if (std::abs(initialCurv.m_kmax) > std::abs(initialCurv.m_kmin))
                {
                    scaleMax = m_featureEnhanceScale;
                }
                if (std::abs(initialCurv.m_kmax) < std::abs(initialCurv.m_kmin))
                {
                    scaleMin = m_featureEnhanceScale;
                }

                m_optimizer.addTargetVertexCurvature(v, scaleMax * initialCurv.m_kmax, scaleMin * initialCurv.m_kmin);
            }
            else if (m_curvatureEditMode == CurvatureEditMode::Clamp)
            {
                m_optimizer.addTargetVertexCurvature(v, std::clamp(initialCurv.m_kmax, m_clampMin, m_clampMax), std::clamp(initialCurv.m_kmin, m_clampMin, m_clampMax));
            }
        }
    }

    // setup boundary curvature targets
    if (m_useBoundaryCurvs)
    {
        VertexProperty<double> savedBoundaryCurvs = m_mesh->get_vertex_property<double>("v:savedBoundaryCurvs");
        for (Vertex v : m_mesh->vertices())
        {
            if (m_mesh->is_boundary(v))
            {
                if (useSavedCurvs)
                {
                    if (savedBoundaryCurvs)
                    {
                        m_optimizer.addTargetBoundaryCurvature(v, savedBoundaryCurvs[v]);
                    }
                }
                else if (m_boundaryCurvatureEditMode == BoundaryCurvatureEditMode::Global)
                {
                    m_optimizer.addTargetBoundaryCurvature(v, m_globalBoundaryCurv);
                }
                else if (m_boundaryCurvatureEditMode == BoundaryCurvatureEditMode::Scale)
                {
                    m_optimizer.addTargetBoundaryCurvature(v, m_scaleBoundaryCurv * m_optimizer.getCurrentBoundaryCurvature(v));
                }
            }
        }
    }

    const auto threadFunc = [&, lockedVertices]()
    {
        m_running = true;

        m_optimizer.compute(lockedVertices);

        m_running = false;
        updateMesh(false);
    };

    std::thread thread(threadFunc);
    thread.detach();
}

void CurvatureOptimizerViewer::drop(int count, const char** paths)
{
    if (count > 0)
    {
        const std::filesystem::path p = paths[0];
        loadMesh(p);
    }
}

void CurvatureOptimizerViewer::draw(const std::string&)
{
    if (m_meshDirty)
    {
        updateMesh(false);
    }

    const CameraState cameraState = { .m_projectionMatrix = projection_matrix_,
                                      .m_modelViewMatrix = modelview_matrix_,
                                      .m_modelViewProjectionMatrix = projection_matrix_ * modelview_matrix_,
                                      .m_screenWidth = width(),
                                      .m_screenHeight = height() };

    m_meshRenderer->draw(cameraState);

    if (m_meshDirty)
    {
        m_meshDirty = false;
        m_meshDirty.notify_all();
    }
}

void CurvatureOptimizerViewer::process_imgui()
{
    BaseViewer::process_imgui();

#ifndef __EMSCRIPTEN__
    if (ImGui::CollapsingHeader("Camera", ImGuiTreeNodeFlags_None))
    {
        if (ImGui::Button("Copy Camera State to Clipboard"))
        {
            const std::string jsonString = getCameraStateAsJson();
            ImGui::SetClipboardText(jsonString.c_str());
        }

        if (ImGui::Button("Apply Camera State from Clipboard"))
        {
            if (ImGui::GetClipboardText())
            {
                const std::string jsonString(ImGui::GetClipboardText());
                setCameraStateAsJson(jsonString);
            }
        }
    }
#endif

    if (!m_mesh->is_empty())
    {
        ImGui::Spacing();

        if (ImGui::CollapsingHeader("Visualization", m_collapseGuiOnStart ? ImGuiTreeNodeFlags_None : ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::PushItemWidth(130);
            float alpha = m_meshRenderer->getMeshAlpha();
            if (ImGui::SliderFloat("Mesh Transparency", &alpha, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_NoRoundToFormat))
            {
                m_meshRenderer->setMeshAlpha(alpha);
            }
            float creaseAngle = m_meshRenderer->getCreaseAngle();
            if (ImGui::SliderFloat("Crease Angle", &creaseAngle, 0.0f, 180.0f, "%.0f"))
            {
                m_meshRenderer->setCreaseAngle(creaseAngle);
            }
            { // draw mode
                const char* current_item = s_drawModeNames[(int) m_drawMode].c_str();
                if (ImGui::BeginCombo("Mesh Draw Mode", current_item))
                {
                    for (size_t i = 0; i < s_drawModeNames.size(); ++i)
                    {
                        const char* item = s_drawModeNames[i].c_str();
                        bool is_selected = (current_item == item);
                        if (ImGui::Selectable(item, is_selected))
                        {
                            m_drawMode = static_cast<DrawMode>(i);
                            updateMeshDrawMode();
                        }
                        if (is_selected)
                        {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
            }
            if (m_drawMode == DrawMode::Curvature)
            { // curvature draw mode
                const char* current_item = s_curvatureDrawModeNames[(int) m_curvatureDrawMode].c_str();
                if (ImGui::BeginCombo("Curvature Draw Mode", current_item))
                {
                    for (size_t i = 0; i < s_curvatureDrawModeNames.size(); ++i)
                    {
                        const char* item = s_curvatureDrawModeNames[i].c_str();
                        bool is_selected = (current_item == item);
                        if (ImGui::Selectable(item, is_selected))
                        {
                            m_curvatureDrawMode = static_cast<CurvatureDrawMode>(i);
                            updateMeshDrawMode();
                        }
                        if (is_selected)
                        {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
            }
            if (ImGui::Checkbox("Show Curvature Directions", &m_showCurvatureDirections))
            {
                updateMeshDrawMode();
            }
            ImGui::PopItemWidth();
        }

        ImGui::Spacing();

        static bool useSavedCurvs = false;
        ImGui::BeginDisabled(m_running);
        if (ImGui::CollapsingHeader("Setup", m_collapseGuiOnStart ? ImGuiTreeNodeFlags_None : ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::PushItemWidth(130);
            {
                static std::array<std::string, 3> optimizerModeNames = { "TargetCurvature", "WillmoreEnergy", "CurvatureVariation" };
                const char* current_item = optimizerModeNames[(int) m_optimizer.getOptimizerMode()].c_str();
                if (ImGui::BeginCombo("Optimizer mode", current_item))
                {
                    for (size_t i = 0; i < optimizerModeNames.size(); ++i)
                    {
                        const char* item = optimizerModeNames[i].c_str();
                        bool is_selected = (current_item == item);
                        if (ImGui::Selectable(item, is_selected))
                        {
                            m_optimizer.setOptimizerMode(static_cast<CurvatureOptimizer::OptimizerMode>(i));
                        }
                        if (is_selected)
                        {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
            }
            double val = m_optimizer.getCurvatureEnergyWeight();
            if (ImGui::InputDouble("Curvature weight", &val, 0.2))
            {
                m_optimizer.setCurvatureEnergyWeight(val);
            }
            val = m_optimizer.getBoundaryCurvatureEnergyWeight();
            if (ImGui::InputDouble("Boundary curvature weight", &val, 0.2))
            {
                m_optimizer.setBoundaryCurvatureEnergyWeight(val);
            }
            val = m_optimizer.getAreaEnergyWeight();
            if (ImGui::InputDouble("Area weight", &val, 0.2))
            {
                m_optimizer.setAreaEnergyWeight(val);
            }
            val = m_optimizer.getConformalEnergyWeight();
            if (ImGui::InputDouble("Conformal weight", &val, 0.2))
            {
                m_optimizer.setConformalEnergyWeight(val);
            }

            ImGui::Spacing();

            ImGui::Checkbox("Lock boundary", &m_lockBoundary);
            if (m_lockBoundary)
            {
                ImGui::Checkbox("Lock two ring", &m_lockTwoRing);
            }

            bool dynamicInitialStepSize = m_optimizer.getDynamicInitialStepSize();
            if (ImGui::Checkbox("Dynamic init. step size", &dynamicInitialStepSize))
            {
                m_optimizer.setDynamicInitialStepSize(dynamicInitialStepSize);
            }

            bool robustDihedralAngles = m_optimizer.getUseRobustDihedralAngles();
            if (ImGui::Checkbox("Robust dihedral angles", &robustDihedralAngles))
            {
                m_optimizer.setUseRobustDihedralAngles(robustDihedralAngles);
            }

            {
                const char* current_item = s_dualAreaModeNames[(int) m_optimizer.getDualAreaMode()].c_str();
                if (ImGui::BeginCombo("Dual area mode", current_item))
                {
                    for (size_t i = 0; i < s_dualAreaModeNames.size(); ++i)
                    {
                        const char* item = s_dualAreaModeNames[i].c_str();
                        bool is_selected = (current_item == item);
                        if (ImGui::Selectable(item, is_selected))
                        {
                            m_optimizer.setDualAreaMode(static_cast<DualAreaMode>(i));
                        }
                        if (is_selected)
                        {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
            }

            {
                const char* current_item = s_vertexCurvatureOperatorNames[(int) m_optimizer.getVertexCurvatureOperator()].c_str();
                if (ImGui::BeginCombo("Vertex curvature operator", current_item))
                {
                    for (size_t i = 0; i < s_vertexCurvatureOperatorNames.size(); ++i)
                    {
                        const char* item = s_vertexCurvatureOperatorNames[i].c_str();
                        bool is_selected = (current_item == item);
                        if (ImGui::Selectable(item, is_selected))
                        {
                            m_optimizer.setVertexCurvatureOperator(static_cast<VertexCurvatureOperator>(i));
                        }
                        if (is_selected)
                        {
                            ImGui::SetItemDefaultFocus();
                        }
                    }
                    ImGui::EndCombo();
                }
            }

            ImGui::PopItemWidth();

            ImGui::Separator();
            ImGui::Spacing();

            if (m_optimizer.getOptimizerMode() == CurvatureOptimizer::OptimizerMode::TargetCurvature)
            {
                if (ImGui::Button("Save current curvs"))
                {
                    saveCurrentCurvatures();
                }
                ImGui::Spacing();
                ImGui::PushItemWidth(130);
                {
                    const char* current_item = s_curvatureEditModeNames[(int) m_curvatureEditMode].c_str();
                    if (ImGui::BeginCombo("Curvature edit mode", current_item))
                    {
                        for (size_t i = 0; i < s_curvatureEditModeNames.size(); ++i)
                        {
                            const char* item = s_curvatureEditModeNames[i].c_str();
                            bool is_selected = (current_item == item);
                            if (ImGui::Selectable(item, is_selected))
                            {
                                m_curvatureEditMode = static_cast<CurvatureEditMode>(i);
                            }
                            if (is_selected)
                            {
                                ImGui::SetItemDefaultFocus();
                            }
                        }
                        ImGui::EndCombo();
                    }
                }
                if (m_curvatureEditMode == CurvatureEditMode::Global)
                {
                    ImGui::InputDouble("Global max curv", &m_globalMaxCurv, 0.2);
                    ImGui::Spacing();
                    ImGui::InputDouble("Global min curv", &m_globalMinCurv, 0.2);
                    ImGui::Spacing();
                }
                else if (m_curvatureEditMode == CurvatureEditMode::Scale)
                {
                    ImGui::InputDouble("Scale max curv", &m_scaleMaxCurv, 0.2);
                    ImGui::Spacing();
                    ImGui::InputDouble("Scale min curv", &m_scaleMinCurv, 0.2);
                    ImGui::Spacing();
                }
                else if (m_curvatureEditMode == CurvatureEditMode::FeatureEnhance)
                {
                    ImGui::InputDouble("Feature upscale", &m_featureEnhanceScale, 0.2);
                }
                else if (m_curvatureEditMode == CurvatureEditMode::Clamp)
                {
                    ImGui::InputDouble("Clamp max", &m_clampMax, 0.2);
                    ImGui::InputDouble("Clamp min", &m_clampMin, 0.2);
                }
                ImGui::Checkbox("Use boundary curv", &m_useBoundaryCurvs);
                if (m_useBoundaryCurvs)
                {
                    {
                        const char* current_item = s_boundaryCurvatureEditModeNames[(int) m_boundaryCurvatureEditMode].c_str();
                        if (ImGui::BeginCombo("Boundary curvature edit mode", current_item))
                        {
                            for (size_t i = 0; i < s_boundaryCurvatureEditModeNames.size(); ++i)
                            {
                                const char* item = s_boundaryCurvatureEditModeNames[i].c_str();
                                bool is_selected = (current_item == item);
                                if (ImGui::Selectable(item, is_selected))
                                {
                                    m_boundaryCurvatureEditMode = static_cast<BoundaryCurvatureEditMode>(i);
                                }
                                if (is_selected)
                                {
                                    ImGui::SetItemDefaultFocus();
                                }
                            }
                            ImGui::EndCombo();
                        }
                    }

                    if (m_boundaryCurvatureEditMode == BoundaryCurvatureEditMode::Global)
                    {
                        ImGui::Spacing();
                        ImGui::InputDouble("Global boundary curv", &m_globalBoundaryCurv, 0.2);
                    }
                    else if (m_boundaryCurvatureEditMode == BoundaryCurvatureEditMode::Scale)
                    {
                        ImGui::Spacing();
                        ImGui::InputDouble("Scale boundary curv", &m_scaleBoundaryCurv, 0.2);
                    }
                }
                ImGui::PopItemWidth();

                ImGui::Separator();
                ImGui::Spacing();

                if (m_mesh->has_vertex_property("v:savedCurvs"))
                {
                    ImGui::Checkbox("Use saved curvs", &useSavedCurvs);
                }
                ImGui::Checkbox("Optimize boundary verts", &m_optimizeBoundary);
            }
            else if (m_optimizer.getOptimizerMode() == CurvatureOptimizer::OptimizerMode::CurvatureVariation)
            {
                ImGui::PushItemWidth(130);
                if (ImGui::InputDouble("pNorm", &m_pNormFairing, 0.2))
                {
                    m_optimizer.setPNormFairing(m_pNormFairing);
                }
                ImGui::PopItemWidth();
            }
        }
        ImGui::Spacing();
        if (ImGui::Button("Optimize"))
        {
            startOptimizer(useSavedCurvs);
        }
        ImGui::EndDisabled();

#ifndef __EMSCRIPTEN__
        ImGui::Separator();
        ImGui::Spacing();
        if (ImGui::Button("Write current mesh"))
        {
            const std::filesystem::path file = "out.ply";
            std::cout << "Write current mesh to: " << file << std::endl;

            std::vector<double> kmax;
            std::vector<double> kmin;
            for (Vertex v : m_mesh->vertices())
            {
                const auto& curv = m_optimizer.getCurrentVertexCurvature(v);
                kmax.push_back(curv.m_kmax);
                kmin.push_back(curv.m_kmin);
            }
            writeMeshPLY(file, *m_mesh, { { "kmax", kmax }, { "kmin", kmin } });
        }
        ImGui::Checkbox("Write mesh sequence", &m_writeSequence);
#endif

        // draw plot window
        if (!m_plots.empty())
        {
            ImGui::Begin("Plots", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoFocusOnAppearing);
            ImGui::SetWindowPos(ImVec2(ImGui::GetMainViewport()->Size.x - ImGui::GetWindowWidth() - 4, 4));
            if (m_collapseGuiOnStart)
            {
                ImGui::SetWindowCollapsed(true, ImGuiCond_Once);
            }

            ImGui::Text("Current energy: %f", m_plots["Energy"].back());
            if (m_optimizer.getOptimizerMode() == CurvatureOptimizer::OptimizerMode::TargetCurvature && m_plots.contains("Curvature Rel Error"))
            {
                ImGui::Text("Current curvature error: %f", m_plots["Curvature Rel Error"].back());
            }

            for (const auto& [name, plot] : m_plots)
            {
                if (ImPlot::BeginPlot(name.c_str(), ImVec2(250, 150), ImPlotFlags_NoLegend))
                {
                    ImPlotStyle& style = ImPlot::GetStyle();
                    style.FitPadding = ImVec2(0.1, 0.1);
                    ImPlot::SetupAxes("", "", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
                    ImPlot::PlotLine(name.c_str(), plot.data(), plot.size());
                    {
                        const double x = plot.size() - 1;
                        const double y = plot.back();
                        const ImVec4 color = ImGui::GetStyle().Colors[ImGuiCol_Button];
                        ImPlot::Annotation(x, y, color, ImVec2(-15, -15), true, "%.4f", plot.back());
                    }
                    ImPlot::EndPlot();
                }
            }

            ImGui::End();
        }
    }
}

void CurvatureOptimizerViewer::mouse(int button, int action, int mods)
{
    if (action == GLFW_PRESS && shift_pressed())
    {
        double x, y;
        cursor_pos(x, y);
        const Vertex v = pickVertex(x, y);
        if (v.is_valid())
        {
            const auto& curv = m_optimizer.getCurrentVertexCurvature(v);
            addLabel(m_mesh->position(v), std::format("{:.4f}, {:.4f}", curv.m_kmax, curv.m_kmin));
        }
    }
    else
    {
        TrackballViewer::mouse(button, action, mods);
    }
}

void CurvatureOptimizerViewer::keyboard(int key, int code, int action, int mods)
{
    if (action != GLFW_PRESS && action != GLFW_REPEAT)
    {
        return;
    }

    switch (key)
    {
    case GLFW_KEY_BACKSPACE: // reload model
    {
        if (!m_meshFile.empty())
        {
            loadMesh(m_meshFile);
        }
        break;
    }
    case GLFW_KEY_SPACE:
    {
        m_drawMode = static_cast<DrawMode>(static_cast<int>(m_drawMode) + 1);
        if (static_cast<int>(m_drawMode) >= static_cast<int>(s_drawModeNames.size()))
        {
            m_drawMode = static_cast<DrawMode>(0);
        }
        updateMeshDrawMode();
        break;
    }
    case GLFW_KEY_N: // toggle normals
    {
        m_meshRenderer->setShowNormals(!m_meshRenderer->getShowNormals());
        break;
    }
    case GLFW_KEY_B: // toggle bounding box
    {
        m_meshRenderer->setShowBoundingBox(!m_meshRenderer->getShowBoundingBox());
        break;
    }
    default:
    {
        TrackballViewer::keyboard(key, code, action, mods);
        break;
    }
    }
}