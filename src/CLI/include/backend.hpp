#ifndef BACKEND_GUARD
#define BACKEND_GUARD

#include <string>
#include <vector>
#include <fstream>
#include <optional>
#include <filesystem>
#include "argparse/argparse.hpp"
#include "argument_exception.hpp"


/**
 * Check to see if a file exists on the system.
 * Technically this will return false if it exists but you don't have permission.
*/
inline bool file_exists(std::string fname) {
    std::fstream f(fname.c_str());
    return f.good();
}

inline std::string absolute_fp(std::string fname) {
    return std::string(std::filesystem::absolute(fname));
}

inline bool can_exec(const std::filesystem::path& fp) {
    auto entry = std::filesystem::directory_entry(fp);
    return (
        entry.exists()
        &&
        std::filesystem::perms::none != (entry.status().permissions() & std::filesystem::perms::owner_exec)
    );
}

static std::vector<std::filesystem::path>& get_paths() {
    std::string path = std::getenv("PATH");
    static std::vector<std::filesystem::path> result;
    if (result.size() > 0) return result;
    result.push_back(std::filesystem::canonical("/proc/self/exe").parent_path());
    result.push_back(std::filesystem::path("."));

    size_t i = 0;
    // Windows uses ";" as a delimeter, but we don't support windows yet.
    while ((i = path.find(":")) != std::string::npos) {
        result.push_back(std::filesystem::path(path.substr(0, i)));
        path.erase(0, i + 1);
    }
    return result;
}

struct FierroBackend {
    std::string name = "";
    std::shared_ptr<argparse::ArgumentParser> command;
    std::optional<std::filesystem::path> exec_path{};
    FierroBackend(std::string name): name(name) {
        exec_path = find();
    }

    std::optional<std::filesystem::path> find() {
        return this->find(get_paths());
    }

    std::optional<std::filesystem::path> find(const std::vector<std::filesystem::path>& paths) {
        for(const auto & path : paths) {
            auto potential_path = path / this->name;
            if (can_exec(potential_path)) return {potential_path};
        }
        return {};
    }

    virtual int invoke() = 0;
};
#endif