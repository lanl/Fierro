#include "backend.hpp"
#include "FierroParallelBackends.hpp"
#include "argparse/argparse.hpp"
#include "argument_exception.hpp"
#include <filesystem>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>

std::vector<std::shared_ptr<FierroBackend>> BACKENDS {
    std::shared_ptr<FierroBackend>(new ParallelExplicit()),
    std::shared_ptr<FierroBackend>(new ParallelImplicit())
};

std::vector<std::shared_ptr<FierroBackend>> find_fierro_backends() {
    auto found = std::vector<std::shared_ptr<FierroBackend>>();
    for(auto backend : BACKENDS) {
        if(backend->exists()) found.push_back(backend);
    }

    return found;
}

/**
 * Execute a function with a particular working directory.
 * 
 * @param dir The working directory to execute the function in.
 *              If the directoy doesn't exist, it will be created.
 *              Will also create parent directories if necessary.
 * @param f   The function to execute.
 * 
 * @return The return value of `f` if there is one. Else void.
*/
template<typename T> 
T with_curdir(std::string dir, std::function<T()> f) {
    T val;
    with_curdir<void>(dir, [&]() -> void { val = f(); });
    return val;
}

/**
 * `with_curdir` template specialization for void functions.
*/
template<>
void with_curdir<void>(std::string dir, std::function<void()> f) {
    auto old_path = std::filesystem::current_path();
    
    dir = std::filesystem::absolute(dir);
    if (!std::filesystem::exists(dir)) std::filesystem::create_directories(dir);

    std::filesystem::current_path(dir);
    f();
    std::filesystem::current_path(old_path);
}

int main(int argc, char** argv) {
    argparse::ArgumentParser parser("fierro");
    parser.add_description("");
    parser.add_argument("-o", "--output")
        .help("output directory for the program")
        .default_value(".");

    auto backends = find_fierro_backends();
    for(auto backend : backends) {
        parser.add_subparser(*backend->command);
    }

    try {
        parser.parse_args(argc, argv);
    } catch(const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }

    for(auto backend : backends) {
        if (parser.is_subcommand_used(*backend->command)) {
            return with_curdir<int>(
                parser.get<std::string>("output"),
                [&]() {
                    try {
                        return backend->invoke();
                    } catch (ArgumentException e) {
                        std::cerr << e.message << std::endl;
                        std::cerr << parser;
                        std::exit(1);
                    }
                }
            );
        }
    }

    // If they didn't select anything, give them some help.
    std::cout << parser << std::endl;
}