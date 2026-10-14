#include "additive_data.hpp"
#include "parse_yaml.hpp"
#include <cstdio>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
const std::string path_file = "laser_input_test.csv";

void write_path(const std::string& text)
{
    std::ofstream file(path_file);
    file << text;
    if (!file) throw std::runtime_error("could not write test path");
}

void expect_failure(const std::function<void()>& operation)
{
    try { operation(); }
    catch (const std::exception&) { return; }
    throw std::runtime_error("invalid laser input was accepted");
}

Laser_t parse(const std::string& text)
{
    Yaml::Node root;
    std::stringstream stream(text);
    Yaml::Parse(root, stream);
    Laser_t laser;
    parse_laser(root, laser);
    return laser;
}

void check(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}

void run_tests()
{
    write_path("# time,x,y,z,power\n\n 0, 1, 2, 3, 10\n2, 3, 6, 7, 30 # end\n");
    auto path = ToolPathInfo::load_laser_path(path_file.c_str());
    mtr::DCArrayKokkos<double> samples(5, 4, "path_samples");
    FOR_ALL(i, 0, 5, {
        const double t = i - 1.0;
        double x;
        double y;
        double z;
        path.get_position(t, x, y, z);
        samples(i, 0) = x;
        samples(i, 1) = y;
        samples(i, 2) = z;
        samples(i, 3) = path.get_power(t);
    }, "sample_laser_path");
    samples.update_host();
    const double expected[5][4] = {{1,2,3,10}, {1,2,3,10}, {2,4,5,20},
                                    {3,6,7,30}, {3,6,7,30}};
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 4; ++j)
            check(samples.host(i, j) == expected[i][j], "device path interpolation mismatch");

    auto absent = parse("dynamic_options:\n  time_final: 1\n");
    check(!absent.configured, "absent laser_options changed defaults");
    auto spherical = parse("laser_options:\n  radius: 0.4\n  model: spherical\n  path: " + path_file + "\n");
    check(spherical.configured && spherical.model == LaserModelType::SPHERICAL &&
          spherical.spherical.radius == 0.4, "spherical configuration mismatch");
    auto goldak = parse("laser_options:\n  path: " + path_file +
        "\n  model: goldak\n  absorptivity: 0.7\n  major_front: 2\n  major_rear: 4\n  minor: 1\n  depth: 0.5\n");
    check(goldak.configured && goldak.model == LaserModelType::GOLDAK &&
          goldak.goldak.major_rear == 4 && goldak.goldak.depth == 0.5,
          "Goldak configuration mismatch");

    const std::string prefix = "laser_options:\n  path: " + path_file + "\n";
    for (const auto& text : {
        std::string("laser_options: invalid\n"),
        prefix + "  model: invalid\n",
        prefix + "  model: spherical\n",
        prefix + "  model: spherical\n  radius: 0\n",
        prefix + "  model: spherical\n  radius: -1\n",
        prefix + "  model: spherical\n  radius: nan\n",
        prefix + "  model: spherical\n  radius: 1junk\n",
        prefix + "  model: spherical\n  radius: 1\n  raduis: 2\n",
        prefix + "  model: goldak\n  absorptivity: 2\n",
        prefix + "  model: goldak\n  absorptivity: 0.5\n  major_front: 1\n",
        std::string("laser_options:\n  model: spherical\n  radius: 1\n")}) {
        expect_failure([&] { parse(text); });
    }

    write_path("0,0,0,0,0\n1,1,0,0,10\n1,2,0,0,20\n2,3,0,0,30\n");
    auto jump = ToolPathInfo::load_laser_path(path_file.c_str());
    check(jump.tool_path_table.linear_interpolation_host(1, ToolPathInfo::x) == 1,
          "duplicate-time transition changed existing exact-time convention");
    check(jump.tool_path_table.linear_interpolation_host(1.5, ToolPathInfo::x) == 2.5,
          "duplicate-time interpolation mismatch");
    write_path("0,1,2,3,0\n");
    auto constant = ToolPathInfo::load_laser_path(path_file.c_str());
    check(constant.tool_path_table.linear_interpolation_host(1, ToolPathInfo::z) == 3,
          "single-row path must remain constant");

    for (const auto& text : {"", "# empty\n", "0,1,2,3\n", "0,1,2,3,4,5\n",
         "0,1,2,3,4,\n", "0,1,,3,4\n", "0,1,2,3,nan\n", "0,1,2,3,inf\n",
         "0,1,2,3,-1\n", "0,1,2,3,4x\n", "1,0,0,0,0\n0,0,0,0,0\n"}) {
        write_path(text);
        expect_failure([&] { ToolPathInfo::load_laser_path(path_file.c_str()); });
    }
    std::remove(path_file.c_str());
    expect_failure([&] { ToolPathInfo::load_laser_path(path_file.c_str()); });
}
} // namespace

int main(int argc, char** argv)
{
    MATAR_INITIALIZE(argc, argv);
    int status = 0;
    try { run_tests(); }
    catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        status = 1;
    }
    std::remove(path_file.c_str());
    MATAR_FINALIZE();
    return status;
}
