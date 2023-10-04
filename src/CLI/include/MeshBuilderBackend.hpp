#ifndef MESH_BUILDER_BACKENDS_H
#define MESH_BUILDER_BACKENDS_H

#include "backend.hpp"
#include "argparse/argparse.hpp"
#include <string>
#include <fstream>

struct MeshBuilderBackend: public FierroBackend {
    MeshBuilderBackend() : FierroBackend("fierro-mesh-builder") {
        this->command = std::shared_ptr<argparse::ArgumentParser>(
            new argparse::ArgumentParser("mesh-builder")
        );
        
        this->command->add_argument("-c", "config")
            .help("The `.yaml` configuration to run fierro with.`") 
            .action(absolute_fp);
        this->command->add_argument("--box-help")
            .help("Set this flag to print out an example box mesh configuration file.")
            .default_value(false)
            .implicit_value(true);
        this->command->add_argument("--cylinder-help")
            .help("Set this flag to print out an example cylinder mesh configuration file.")
            .default_value(false)
            .implicit_value(true);;
        this->command->add_description("Build rectangular or cylindrical mesh in 2 or 3 dimensions. Useful for setting up test problems.");
    }

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
        if (this->command->is_used("box-help") || this->command->is_used("cylinder-help"))
            return {};

        if (!parser->is_used("config"))
            msg = "A value for config is required";
        else if (!file_exists(parser->get<std::string>("config"))) {
            msg = "Unable to find configuration: " + parser->get<std::string>("config");
        }
        
        if (msg.length() > 0) {
            return std::optional(msg);
        }
        return {};
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

        std::string sys_command;
        if (this->command->get<bool>("box-help")) 
            sys_command = this->exec_path.value().string() + " Box";
        else if (this->command->get<bool>("cylinder-help"))
            sys_command = this->exec_path.value().string() + " Cylinder";
        else {
            std::string input_file = this->command->get<std::string>("config");
            sys_command = this->exec_path.value().string() + " " + input_file;
        }

        return system(sys_command.c_str());
    }

    void add_common_options() {
        this->command->add_argument("config")
            .help("The `.yaml` configuration to run fierro with.`")
            .action(absolute_fp);
    }
};

#endif