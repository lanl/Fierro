#ifndef FIERRO_PARALLEL_BACKENDS
#define FIERRO_PARALLEL_BACKENDS
#include "backend.hpp"
#include "argparse/argparse.hpp"
#include <string>
#include <fstream>

/**
 * Check to see if a file exists on the system.
 * Technically this will return false if it exists but you don't have permission.
*/
bool file_exists(std::string fname) {
    std::fstream f(fname.c_str());
    return f.good();
}

std::string absolute_fp(std::string fname) {
    return std::string(std::filesystem::absolute(fname));
}

struct ParallelBackend: public FierroBackend {
    ParallelBackend(std::string name) 
        : FierroBackend(name) { }

    /** 
     * Check that the arguments of the CLI are valid and make sense.
     * 
     * Checks for explicit/implicit logical consistency.
     * Checks that the mesh file is real and readable.
     *      Does not check that it is correctly formatted.
     * 
     * @return An optional error. If the empty, the arguments are valid.
     * 
    */
    std::optional<std::string> validate_arguments() {
        std::string msg = "";
        auto parser = this->command;

        if (parser->is_used("mesh_file") && !file_exists(parser->get<std::string>("mesh_file"))) {
            msg = "Unable to find mesh file: " + parser->get<std::string>("mesh_file");
        } else if (parser->is_used("config") && !file_exists(parser->get<std::string>("config"))) {
            msg = "Unable to find configuration: " + parser->get<std::string>("config");
        } else if (parser->is_used("config") == parser->is_used("mesh_file")) {
            msg = "Use exactly one of `--config` or `--mesh_file.`";
        }
        
        if (msg.length() > 0) {
            return std::optional(msg);
        }
        return std::nullopt;
    }
    
    /**
     * Invoke the backend executable.
     * Will throw an ArgumentException if the arguments aren't valid.
     * 
     * @return Return code of the executable.
    */
    int invoke() {
        auto err = this->validate_arguments();
        if (err.has_value()) throw new ArgumentException(err.value());

        std::string input_file;
        if (this->command->is_used("config")) {
            input_file = this->command->get<std::string>("config");
        } else {
            input_file = this->command->get("mesh_file");
        }
        
        std::string sys_command = this->name + " " + input_file;

        // If the user gives us a number of processes,
        // we should set up a couple of environment variables and invoke mpi.
        if (this->command->is_used("np")) {
            system("export OMP_PROC_BIND=spread");
            system("export OMP_NUM_THREADS=1   ");
            system("export OMP_PLACES=threads  ");
            sys_command = "mpirun -np " + this->command->get<std::string>("np") + " --bind-to core " + sys_command;
        }

        return system(sys_command.c_str());
    }

    void add_common_options() {
        this->command->add_argument("-m", "--mesh_file")
            .help("The `.geo` file to run fierro with. Mutually exclusive with `--config.`")
            .action(absolute_fp);

        this->command->add_argument("-c", "--config")
            .help("The `.yaml` configuration to run fierro with. Mutually exclusive with `--mesh_file.`")
            .action(absolute_fp);

        this->command->add_argument("-np")
            .help("Number of processes to run. If set, we will invoke `mpirun -np {v}`")
            .default_value("1"); // This is an int, but were just going to put it back into a str anyway.
    }
};

struct ParallelExplicit: public ParallelBackend {
    ParallelExplicit() : ParallelBackend("fierro-parallel-explicit") {
        this->command = std::shared_ptr<argparse::ArgumentParser>(
            new argparse::ArgumentParser("parallel-explicit")
        );
        this->add_common_options();
        this->command->add_description("");
    }
};

struct ParallelImplicit: public ParallelBackend {
    ParallelImplicit() : ParallelBackend("fierro-parallel-implicit") {
        this->command = std::shared_ptr<argparse::ArgumentParser>(
            new argparse::ArgumentParser("parallel-implicit"));
        this->add_common_options();
        this->command->add_description("");
    }
};
#endif