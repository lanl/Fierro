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
std::cout << "Required arguments:\n"
             "  -x, --x-dim    <value>    x dimension of RVE\n"
             "  -y, --y-dim    <value>    y dimension of RVE\n"
             "  -z, --z-dim    <value>    z dimension of RVE\n"
             "  -f, --infile   <path>     input file (full path)\n"
             "\n"
             "Optional arguments:\n"
             "  -m, --micro-filetype <value>    microstructure file type (0: ASCII, 1: HDF5)\n"
             "  -e, --euler-angles   <path>     full path to euler angles dataset in HDF5 file\n"
             "  -g, --feature-ids    <path>     full path to feature ids dataset in HDF5 file\n"
             "  -p, --phase-ids      <path>     full path to phase ids dataset in HDF5 file\n"
             "  -u, --euler-unit     <value>    unit of euler angles in HDF5 file (0: degree, 1: radian)\n"
             "  -h, --help                      print this message\n"
             "\n"
            "Examples:\n"
             "  Los Alamos FFT ASCII microstructure filetype:\n"
             "    ./executable -x 16 -y 16 -z 16 -f <inputfile_path>\n"
             "    ./executable --x-dim=16 --y-dim=16 --z-dim=16 --infile=<inputfile_path>\n"
             "  HDF5 microstructure filetype:\n"
             "    ./executable -x 16 -y 16 -z 16 -f <inputfile_path> -m 1 -e <eulerAngles_path> -g <featureIds_path> -p <phaseIds_path> -u 1\n";
    
   exit(1);
}
