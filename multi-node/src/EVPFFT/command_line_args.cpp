#include "command_line_args.h"
#include <getopt.h>
#include <iostream>

CommandLineArgs::CommandLineArgs()
{
}

void CommandLineArgs::parse_command_line(int argc, char *argv[])
{
    const char* const short_opts = "x:y:z:f:m:e:g:p:u:h";
    const option long_opts[] = {{"x-dim", required_argument, nullptr, 'x'},
                                {"y-dim", required_argument, nullptr, 'y'},
                                {"z-dim", required_argument, nullptr, 'z'},
                                {"infile", required_argument, nullptr, 'f'},
                                {"micro-filetype", required_argument, nullptr, 'm'},
                                {"euler-angles", required_argument, nullptr, 'e'},
                                {"feature-ids", required_argument, nullptr, 'g'},
                                {"phase-ids", required_argument, nullptr, 'p'},
                                {"euler-unit", required_argument, nullptr, 'u'},
                                {"help", no_argument, nullptr, 'h'},
                                {nullptr, no_argument, nullptr, 0}};

    while (true)
    {
        const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

        if (-1 == opt)
            break;

        switch (opt)
        {
        case 'x': 
          nn[0] = atoi(optarg);
          break;

        case 'y': 
          nn[1] = atoi(optarg); 
          break;

        case 'z': 
          nn[2] = atoi(optarg);
          break;

        case 'f':
          input_filename = optarg;
          break;

        case 'm':
          micro_filetype = atoi(optarg);
          break;

        case 'e':
          EulerAnglesFullPath = optarg;
          break;

        case 'g':
          FeatureIdsFullPath = optarg;
          break;

        case 'p':
          PhaseIdsFullPath = optarg;
          break;

        case 'u':
          EulerAnglesUnit = atoi(optarg);
          break;

        case 'h': // -h or --help
          print_help(); 
          break;

        case '?': // Unrecognized option
          print_help();
          break;

        default:
          print_help();
          break;
        }  // end switch
    }  // end while

    check_cmd_args();
}

void CommandLineArgs::check_cmd_args()
{
    if (nn[0] == 0 or 
        nn[1] == 0 or 
        nn[2] == 0 or
        micro_filetype == -1 or 
        input_filename.size() == 0)
    {
      print_help();
    }

    if (micro_filetype == 1) {
      if (EulerAnglesFullPath.size() == 0 or
          FeatureIdsFullPath.size() == 0 or
          PhaseIdsFullPath.size() == 0 or
          EulerAnglesUnit == -1)
      {
        print_help();
      }
    }
}

void CommandLineArgs::print_help()
{    
    std::cout << "  Required arguments are:\n"
                 "    long          short  arg      description\n"
                 "  --x-dim          -x    1    x dimension of RVE\n"
                 "  --y-dim          -y    1    y dimension of RVE\n"
                 "  --z-dim          -z    1    z dimension of RVE\n"
                 "  --infile         -f    1    input file (full path)\n"
                 "  --micro-filetype -m    1    microstructure file type (ASCII:0, HDF5:1)\n"
                 "  Optional arguments are:\n"
                 "  --euler-angles   -e    1    full path to euler angles dataset in hdf5 file\n"
                 "  --feature-ids    -g    1    full path to feature ids dataset in hdf5 file\n"
                 "  --phase-ids      -p    1    full path to phase ids dataset in hdf5 file\n"
                 "  --euler-unit     -u    1    unit of euler angles in hdf5 file (degree:0, radian:1)\n"
                 "  --help           -h    0    print this message\n"
                 "  More info:\n"
                 "    Note: options with arg=1 require values.\n"
                 "    Example for Los Alamos FFT ASCII microstructure filetype:\n"
                 "        ./executable -x 16 -y 16 -z 16 -f <inputfile_path> -m 0\n"
                 "        ./executable --x-dim=16 --y-dim=16 --z-dim=16 --infile=<inputfile_path> --micro-filetype=0\n"
                 "    Example for hdf5 microstructure filetype:\n"
                 "        ./executable -x 16 -y 16 -z 16 -f <inputfile_path> -m 1 -e <eulerAngles_path> -g <featureIds_path> -p <phaseIds_path> -u 1\n";
    exit(1);
}
