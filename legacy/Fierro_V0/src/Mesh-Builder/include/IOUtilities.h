#pragma once
#include <filesystem>
#include <string>

namespace IOUtilities {
    /**
     * Ensures that a directory exists. 
    */
    inline void mkdir(std::string path) {
        path = std::filesystem::absolute(path);
        if (!std::filesystem::exists(path)) std::filesystem::create_directories(path);
    }
}