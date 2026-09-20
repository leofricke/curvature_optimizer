#pragma once

#include <pmp/types.h>

using namespace pmp;

Color hsvColor(float h, float s, float v);
std::vector<Color> randomColors(size_t n, int seed = -1);