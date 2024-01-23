# Fiero Command Line Interface
This CLI package is a Fierro command line interface that is intended to act as a consolidated front-end for our Fierro backend tools. It provides a consistent interface and framework for exposing tools to the user.

## Front-End
The front end here is built off of [argparse](https://github.com/p-ranav/argparse), which implements a language agnostic standard for simple command line interface design. It defines a syntax for programs from the command line with suitable helptext/defaults/post-processing.

The Fierro-CLI implements two things: 
1. Output directory -- A directory to run the program in. Outputs of the program will by default be put here if some other backend specific option is not specified.
2. Subcommands -- a list of argparse "subcommands" that each represent an available backend tool.

Example:
```
$ fierro -h

Usage: fierro [--help] [--version] [--output VAR] {mesh-builder,parallel-explicit,parallel-implicit,voxelizer}

Optional arguments:
  -h, --help            shows help message and exits
  -v, --version         prints version information and exits
  -o, --output          output directory for the program [default: "."]

Subcommands:
  mesh-builder    Build rectangular or cylindrical mesh in 2 or 3 dimensions. Useful for setting up test problems.
  parallel-explicit Use an explicit solution scheme to step the system.
  parallel-implicit Use an implicit solution scheme to step the system.
  voxelizer       Take a mesh defined by an STL and turn it into a dense voxel grid represented by a VTK structured mesh. The extent of the grid is selected by taking the minimum and maximum extent of the mesh along each dimenison. Note: Only supports binary STL mesh inputs
```

## Backends
This Fierro-CLI package is designed to operate decoupled from backends. Each registered backend is an executable that may or may not be present in the system.

### Registering Backends
To add a new backend, you only need to do a couple of steps.
#### 1. Implementing NewBackend.hpp
You should start by copying one of the existing backend.hpp files and implementing your new backend. We will use `VoxelizerBackend.hpp` as an example.

You should derive from `FierroBackend` and tell the base class which executable name to look for
```c++
#include "backend.hpp"
struct VoxelizerBackend: public FierroBackend {
    VoxelizerBackend() : FierroBackend("fierro-voxelizer") {
        ...
    }
    ...
};
```

Then, in the constructor, setup your argparse command

```c++
    
struct VoxelizerBackend: public FierroBackend {
    VoxelizerBackend() : FierroBackend("fierro-voxelizer") {

        this->command = std::shared_ptr<argparse::ArgumentParser>(
            new argparse::ArgumentParser("voxelizer") // <--- This is the name that shows up under "Subcommands:" in the CLI
        );
        
        // Here you can add your arguments with their helptext.
        this->command->add_argument("input-file")
            .help("The binary STL file to voxelize.") 
            .action(absolute_fp); // <--  You can also add a post-processing command. In this case we convert the input filepath to an absolute one.
        this->command->add_argument("output-name")
            .help("The name of the VTK structured mesh output file.");
        this->command->add_argument("x")
            .help("The number of voxels along the X direction (int).");
        this->command->add_argument("y")
            .help("The number of voxels along the Y direction (int).");
        this->command->add_argument("z")
            .help("The number of voxels along the Z direction (int).");
        this->command->add_argument("use-index-space")
            .help("Wether or not to output the VTK in index-space coordinates.");

        // Add a description for your subcommand
        this->command->add_description(
            "Take a mesh defined by an STL and turn it into a dense voxel grid represented by a VTK structured mesh. "
            "The extent of the grid is selected by taking the minimum and maximum extent of the mesh along each dimenison. "
            "Note: Only supports binary STL mesh inputs."
        );
    }
    ...
};
```

Lastly you tell the CLI framework how to actually invoke the backend command. This is mostly just building a system command string and invoking it, but I suppose you could do whatever you want.
```c++

struct VoxelizerBackend: public FierroBackend {
    ...
    
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
```

You probably also want to add some validation to your input arguments. The voxelizer does this with a serparate function:

```c++
struct VoxelizerBackend: public FierroBackend {
    ...
    
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

        // Here we just check if the input file is an actual file.
        if (!file_exists(parser->get<std::string>("input-file"))) {
            msg = "Unable to find input file: " + parser->get<std::string>("input-file");
        }
        
        if (msg.length() > 0) {
            return std::optional(msg);
        }
        return {};
    }
    ...
};
```

Now you have a well defined backend.


#### 2. Plug it into the CLI framework
Head over to main.cpp and add your backend to the list of backends:
```c++
std::vector<std::shared_ptr<FierroBackend>> BACKENDS {
    std::shared_ptr<ParallelExplicit>(),
    std::shared_ptr<ParallelImplicit>(),
    std::shared_ptr<MeshBuilderBackend>(),
    std::shared_ptr<VoxelizerBackend>(),
    //std::shared_ptr<YourBackend>() <-- Add yours here
};

```

Now you have everything set up, and it will automatically be picked up if its compiled and in the system.