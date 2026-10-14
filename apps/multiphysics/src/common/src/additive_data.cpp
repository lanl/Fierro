#include "additive_data.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// File I/O and host path sampling stay out of device-facing headers.
ToolPathInfo ToolPathInfo::load_laser_path(const char* path_name) {
    const std::string filename(path_name);
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open laser path: " + filename);
    struct Row { double values[5]; };
    std::vector<Row> rows;
    std::string line;
    size_t line_number = 0;
    while (std::getline(file, line)) {
        ++line_number;
        line = line.substr(0, line.find('#'));
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        const std::string location = filename + ":" + std::to_string(line_number);
        Row entry{};
        double* row = entry.values;
        std::istringstream stream(line);
        std::string field;
        size_t column = 0;
        while (std::getline(stream, field, ',')) {
            if (column == 5)
                throw std::runtime_error(location + ": expected five CSV fields");
            const auto end = field.find_last_not_of(" \t\r\n");
            if (end == std::string::npos)
                throw std::runtime_error(location + ": empty CSV field");
            field.erase(end + 1);
            size_t consumed = 0;
            double value;
            try { value = std::stod(field, &consumed); }
            catch (const std::exception&) {
                throw std::runtime_error(location + ": invalid numeric CSV field");
            }
            if (consumed != field.size() || !Kokkos::isfinite(value))
                throw std::runtime_error(location + ": expected a finite number");
            row[column++] = value;
        }
        const auto last = line.find_last_not_of(" \t\r\n");
        if (column != 5 || line[last] == ',')
            throw std::runtime_error(location + ": expected five CSV fields");
        if (!rows.empty() && row[0] < rows.back().values[0])
            throw std::runtime_error(location + ": times must be nondecreasing");
        if (row[4] < 0.0)
            throw std::runtime_error(location + ": laser power must be nonnegative");
        rows.push_back(entry);
    }
    if (file.bad()) throw std::runtime_error("Failed reading laser path: " + filename);
    if (rows.empty()) throw std::runtime_error("Laser path is empty: " + filename);
    ToolPathInfo path(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto& r = rows[i].values;
        path.set_data_point(i, r[0], r[1], r[2], r[3], r[4]);
    }
    path.update_device();
    return path;
}

LaserSource_t sample_laser_source(const Laser_t& laser, double time)
{
    const auto& table = laser.tool_path_info.tool_path_table;
    if (!laser.configured || table.num_rows == 0)
        throw std::runtime_error("Cannot sample an unconfigured laser");
    LaserSource_t source;
    source.model = laser.model;
    source.spherical = laser.spherical;
    source.goldak = laser.goldak;
    source.x = table.linear_interpolation_host(time, ToolPathInfo::x);
    source.y = table.linear_interpolation_host(time, ToolPathInfo::y);
    source.z = table.linear_interpolation_host(time, ToolPathInfo::z);
    source.power = table.linear_interpolation_host(time, ToolPathInfo::power);

    // Use the actual piecewise-linear path tangent, independent of timestep.
    // At a waypoint use the incoming segment, as for table interpolation.
    size_t segment = table.num_rows;
    for (size_t i = 0; i + 1 < table.num_rows; ++i) {
        if (table.data.host(i + 1, ToolPathInfo::time) <=
            table.data.host(i, ToolPathInfo::time)) continue;
        segment = i;
        if (time <= table.data.host(i + 1, ToolPathInfo::time)) break;
    }
    if (segment < table.num_rows) {
        const double dx = table.data.host(segment + 1, ToolPathInfo::x) -
                          table.data.host(segment, ToolPathInfo::x);
        const double dy = table.data.host(segment + 1, ToolPathInfo::y) -
                          table.data.host(segment, ToolPathInfo::y);
        const double distance = Kokkos::hypot(dx, dy);
        if (distance > 0.0) {
            source.direction_x = dx / distance;
            source.direction_y = dy / distance;
        }
    }
    return source;
}
