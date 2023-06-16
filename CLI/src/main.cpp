#include "argparse/argparse.hpp"
#include <fstream>
#include <stdlib.h>

bool file_exists(std::string fname) {
    std::fstream f(fname.c_str());
    return f.good();
}

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

    std::string mesh_file = parser.get("mesh_file");
    std::string command = "";
    
    if (parser.get<bool>("explicit")) {
        command = "fierro-parallel-explicit " + mesh_file;
    }
    if (parser.get<bool>("implicit")) {
        command = "fierro-parallel-implicit " + mesh_file;
    }

    // If the user gives us a number of processes,
    // we should set some things and invoke mpi.
    if (parser.is_used("np")) {
        command = "mpirun -np " + parser.get<std::string>("np") + " --bind-to core " + command;
        system("export OMP_PROC_BIND=spread");
        system("export OMP_NUM_THREADS=1");
        system("export OMP_PLACES=threads");
    }

    return system(command.c_str());    
}