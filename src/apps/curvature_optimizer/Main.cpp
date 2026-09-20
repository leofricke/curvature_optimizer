#include "CurvatureOptimizerViewer.h"

#include <args.hxx>

int main(int argc, char** argv)
{
    args::ArgumentParser parser("Curvature Optimizer.", "");
    args::HelpFlag help(parser, "help", "Display this help menu", { 'h', "help" });

    args::Positional<std::string> meshFile(parser, "mesh", "The mesh file to load");
    args::Flag hideGui(parser, "Hide Gui", "Hide Gui on start", { "g" });
    args::Flag collapseGui(parser, "Collapse Gui", "Collapse Gui on start", { "c" });
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Help&)
    {
        std::cout << parser;
        return 0;
    }
    catch (const args::ParseError& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    CurvatureOptimizerViewer viewer("Curvature Optimizer", 1280, 720, hideGui.Get(), collapseGui.Get());
    if (meshFile)
    {
        viewer.loadMesh(meshFile.Get());
    }

    return viewer.run();
}