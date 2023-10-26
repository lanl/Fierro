#ifndef FIERRO_PARALLEL_BACKENDS
#define FIERRO_PARALLEL_BACKENDS
#include "backend.hpp"
#include "argparse/argparse.hpp"
#include <string>
#include <fstream>

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

        if (!file_exists(parser->get<std::string>("config"))) {
            msg = "Unable to find configuration: " + parser->get<std::string>("config");
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
        if (!exec_path.has_value())
            throw std::runtime_error("Cannot invoke executable without absolute path.");
        auto err = this->validate_arguments();
        if (err.has_value()) throw ArgumentException(err.value());

        std::string input_file = this->command->get<std::string>("config");
        std::string sys_command = this->exec_path.value().string() + " " + input_file;

        return system(sys_command.c_str());
    }

    void add_common_options() {
        this->command->add_argument("config")
            .help("The `.yaml` configuration to run fierro with.`")
            .action(absolute_fp);
    }
};

struct ParallelExplicit: public ParallelBackend {
    ParallelExplicit() : ParallelBackend("fierro-parallel-explicit") {
        this->command = std::shared_ptr<argparse::ArgumentParser>(
            new argparse::ArgumentParser("parallel-explicit")
        );
        this->add_common_options();
        this->command->add_description("Use an explicit solution scheme to step the system.");
    }
};

struct ParallelImplicit: public ParallelBackend {
    ParallelImplicit() : ParallelBackend("fierro-parallel-implicit") {
        this->command = std::shared_ptr<argparse::ArgumentParser>(
            new argparse::ArgumentParser("parallel-implicit")
        );
        this->add_common_options();
        this->command->add_description("Use an implicit solution scheme to step the system.");
    }
};
#endif