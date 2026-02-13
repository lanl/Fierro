#pragma once

#include <array>
#include <string>

struct CommandLineArgs
{
    std::string input_filename; // name and path of input file (./fft.in)

    // micro_filetype = 0 means classic Los Alamos FFT ASCII microstructure filetype
    // micro_filetype = 1 means hdf5 microstructure filetype
    int micro_filetype = 0; // default to ASCII

    // if micro_filetype = 1 then information about the full path to locate 
    // the euler angles, featureId, and phases in the hdf5 must be specified;
    // also the euler angles units (degree:0 or radians:1) must be specified
    std::string EulerAnglesFullPath;
    std::string FeatureIdsFullPath;
    std::string PhaseIdsFullPath;
    int EulerAnglesUnit = -1; // set to -1 initially to check that user specified it

    CommandLineArgs();
    void parse_command_line(int argc, char *argv[]);
    void check_cmd_args();
    void print_help();
};
