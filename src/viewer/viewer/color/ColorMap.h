#pragma once

#include <pmp/types.h>

class ColorMap
{
public:
    ColorMap();
    ColorMap(const std::vector<pmp::Color>& values);

    pmp::Color value(float val) const;

private:
    const std::vector<pmp::Color> m_values;
};

class ViridisColorMap : public ColorMap
{
public:
    ViridisColorMap();
};