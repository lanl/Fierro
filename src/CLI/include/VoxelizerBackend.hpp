#ifndef VOXELIZER_BACKENDS_H
#define VOXELIZER_BACKENDS_H

#include "backend.hpp"
#include "argparse/argparse.hpp"
#include <string>
#include <fstream>

struct VoxelizerBackend: public FierroBackend {
    VoxelizerBackend() : FierroBackend("fierro-voxelizer") {
        this->command = std::shared_ptr<argparse::ArgumentParser>(
            new argparse::ArgumentParser("voxelizer")
        );
        
        this->command->add_argument("input-file")
            .help("The binary STL file to voxelize.") 
            .action(absolute_fp);
        this->command->add_argument("output-name")
            .help("The name of the VTK structured mesh output file.");
        this->command->add_argument("x")
            .help("The number of voxels along the X direction (int).");
        this->command->add_argument("y")
            .help("The number of voxels along the Y direction (int).");
        this->command->add_argument("z")
            .help("The number of voxels along the Z direction (int).");
        this->command->add_description(
            "Take a mesh defined by an STL and turn it into a dense voxel grid represented by a VTK structured mesh. "
            "The extent of the grid is selected by taking the minimum and maximum extent of the mesh along each dimenison. "
            "Note: Only supports binary STL mesh inputs."
        );
        this->command->add_argument("use-index-space")
            .help("Wether or not to output the VTK in index-space coordinates.");
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

        if (!file_exists(parser->get<std::string>("input-file"))) {
            msg = "Unable to find input file: " + parser->get<std::string>("input-file");
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
        sys_command = this->exec_path.value().string()
                        + " " + this->command->get("input-file")
                        + " " + this->command->get("output-name")
                        + " " + this->command->get("x")
                        + " " + this->command->get("y")
                        + " " + this->command->get("z")
                        + " " + this->command->get("use-index-space");

        return system(sys_command.c_str());
    }
};

#endif