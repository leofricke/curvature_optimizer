#include "ColorUtil.h"

#include <vector>
#include <random>

Color hsvColor(float h, float s, float v)
{
    const float fC = v * s; // Chroma
    const float hPrime = std::fmod(h / 60.0, 6);
    const float fX = fC * (1 - std::fabs(std::fmod(hPrime, 2) - 1));
    const float fM = v - fC;

    Color color;
    if (0 <= hPrime && hPrime < 1)
        color = Color(fC, fX, 0);
    else if (1 <= hPrime && hPrime < 2)
        color = Color(fX, fC, 0);
    else if (2 <= hPrime && hPrime < 3)
        color = Color(0, fC, fX);
    else if (3 <= hPrime && hPrime < 4)
        color = Color(0, fX, fC);
    else if (4 <= hPrime && hPrime < 5)
        color = Color(fX, 0, fC);
    else if (5 <= hPrime && hPrime < 6)
        color = Color(fC, 0, fX);
    else
        color = Color(0, 0, 0);

    color += Color(fM);
    return color;
}

std::vector<Color> randomColors(size_t n, int seed)
{
    std::vector<Color> colors;
    colors.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        const float c = i / static_cast<double>(n - 1);
        colors.push_back(hsvColor(255.0f * c, 1.0f, 1.0f));
    }

    std::random_device rd;
    std::default_random_engine rng(rd());
    if (seed >= 0)
    {
        rng.seed(seed);
    }
    std::shuffle(colors.begin(), colors.end(), rng);

    return colors;
}