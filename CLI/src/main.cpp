#include "argparse/argparse.hpp"
#include <filesystem>
#include <fstream>
#include <stdlib.h>

/**
 * Check to see if a file exists on the system.
 * Technically this will return false if it exists but you don't have permission.
*/
bool file_exists(std::string fname) {
    std::fstream f(fname.c_str());
    return f.good();
}

/**
 * Check that the arguments of the CLI are valid and make sense.
 * 
 * Checks for explicit/implicit logical consistency.
 * Checks that the mesh file is real and readable.
 *      Does not check that it is correctly formatted.
 * 
 * @param parser ArgumentParser *after* you have called `parser.parse_args(...).`
 * @return An optional error. If the empty, the arguments are valid.
 * 
*/
std::optional<std::string> validate_arguments(const argparse::ArgumentParser& parser) {
    std::string msg = "";

    bool use_explicit = parser.get<bool>("explicit");
    bool use_implicit = parser.get<bool>("implicit");

    if (use_explicit && use_implicit) {
        msg = "Only set one of `--explicit` or `--implicit`";
    } else if (!(use_explicit || use_implicit)) {
        msg = "Requires one of `--explicit` or `--implicit` to be set";
    } else if (!file_exists(parser.get<std::string>("mesh_file"))) {
        msg = "Unable to find mesh file: " + parser.get<std::string>("mesh_file");
    }

    if (msg.length() > 0) {
        return std::optional(msg);
    }
    return std::nullopt;
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

    parser.add_argument("mesh_file")
        .help("The `.geo` file to run fierro on.")
        .required();

    parser.add_argument("-e", "--explicit")
        .help("Use an explicit solution scheme to step the system.")
        .default_value(false)
        .implicit_value(true);
        
    parser.add_argument("-i", "--implicit")
        .help("Use an implicit solution scheme to step the system.")
        .default_value(false)
        .implicit_value(true);
    
    parser.add_argument("-np")
        .help("Number of processes to run. If set, we will invoke `mpirun -np {v}`")
        .default_value("1"); // This is an int, but were just going to put it back into a str anyway.

    parser.add_argument("-o")
        .help("Output folder. ")
        .default_value(".");

    try {
        parser.parse_args(argc, argv);
    } catch(const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        std::exit(1);
    }

    auto err = validate_arguments(parser);
    if (err.has_value()) {
        std::cerr << "Argument Error: " << err.value() << std::endl;
        std::exit(1);
    }

    // We might cd later, so make this an absolute file path first.
    std::string mesh_file = std::filesystem::absolute(parser.get("mesh_file"));
    std::string command = "";
    
    if (parser.get<bool>("explicit")) {
        command = "fierro-parallel-explicit " + mesh_file;
    }
    if (parser.get<bool>("implicit")) {
        command = "fierro-parallel-implicit " + mesh_file;
    }

    // If the user gives us a number of processes,
    // we should set up a couple of environment variables and invoke mpi.
    if (parser.is_used("np")) {
        system("export OMP_PROC_BIND=spread");
        system("export OMP_NUM_THREADS=1   ");
        system("export OMP_PLACES=threads  ");
        command = "mpirun -np " + parser.get<std::string>("np") + " --bind-to core " + command;
    }

    return with_curdir<int>(
        parser.get<std::string>("o"),
        [&]() {
            return system(command.c_str());
        }
    );
}