#pragma once

#include <filesystem>

static std::filesystem::path getTestFile(const std::string& file)
{
    std::filesystem::path path(TEST_FILES_DIR);
    path /= file;
    return path;
}