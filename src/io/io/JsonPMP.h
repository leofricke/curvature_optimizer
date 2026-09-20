#pragma once

#include <pmp/mat_vec.h>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

// pmp vector/matrix serialization support
namespace pmp
{

template <typename T, int M, int N>
void to_json(json& j, const Matrix<T, M, N>& m)
{
    j = json::array();

    for (size_t r = 0; r < M; ++r)
    {
        for (size_t c = 0; c < N; ++c)
        {
            j.push_back(m(r, c));
        }
    }
}

template <typename T, int M, int N>
void from_json(const json& j, Matrix<T, M, N>& m)
{
    for (size_t r = 0; r < M; ++r)
    {
        for (size_t c = 0; c < N; ++c)
        {
            m(r, c) = j[r * N + c];
        }
    }
}

} // namespace pmp
