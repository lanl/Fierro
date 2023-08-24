#ifndef BACKEND_GUARD
#define BACKEND_GUARD

#include <string>
#include <vector>
#include <filesystem>
#include "argparse/argparse.hpp"
#include "argument_exception.hpp"

bool can_exec(const std::filesystem::path& fp) {
    auto entry = std::filesystem::directory_entry(fp);
    return (
        entry.exists()
        &&
        std::filesystem::perms::none != (entry.status().permissions() & std::filesystem::perms::owner_exec)
    );
}

struct FierroBackend {
    std::string name = "";
    std::shared_ptr<argparse::ArgumentParser> command;
    FierroBackend(std::string name): name(name) { }

    bool exists() {
        std::string path = std::getenv("PATH");
        auto result = std::vector<std::filesystem::path>();
        result.push_back(std::filesystem::canonical("/proc/self/exe").parent_path());
        result.push_back(std::filesystem::path("."));

        size_t i = 0;
        // Windows uses ";" as a delimeter, but we don't support windows yet.
        while ((i = path.find(":")) != std::string::npos) {
            result.push_back(std::filesystem::path(path.substr(0, i)));
            path.erase(0, i + 1);
        }

        return this->exists(result);
    }

    bool exists(const std::vector<std::filesystem::path>& paths) {
        for(const auto & path : paths) {
            auto potential_path = path / this->name;
            if (can_exec(potential_path)) return true;
        }
        return false;
    }

    virtual int invoke() = 0;
};
#endif