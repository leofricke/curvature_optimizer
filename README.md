# Robust Curvature-Driven Mesh Optimization

This repository contains the source code of the paper "Robust Curvature-Driven Mesh Optimization" by Leonard Fricke, Ulrich Schwanecke, and Mario Botsch.

## Requirements

- CMake
- C++20 compiler (tested with GCC, Clang, and MSVC)
- Intel MKL (optional, for PARDISO solver)

## Building

    mkdir build && cd build && cmake .. && make -j

- Dependencies are downloaded automatically by cmake using Fetch-Content
- To build with Intel MKL for the PARDISO solver set `MKL_DIR=/opt/intel/oneapi/mkl/latest`

### CMake options

| Option | Default | Description |
| --- | --- | --- |
| `BUILD_HEADLESS` | `OFF` | Build without the GUI |
| `BUILD_TESTS` | `ON` | Build the test suite |

## Usage

Run the viewer with a mesh file:

    ./build/bin/curvature_optimizer path/to/mesh.obj

## Testing

    ./build/bin/test_runner

## Citation

    @inproceedings{fricke_2026,
        author = {Fricke, Leonard and Schwanecke, Ulrich and Botsch, Mario},
        title = {{Robust Curvature-Driven Mesh Optimization}},
        booktitle = {Vision, Modeling, and Visualization},
        editor = {Thomas Leimkühler, Tino Weinkauf and Marcel Campen},
        year = {2026},
        publisher = {The Eurographics Association}
    }