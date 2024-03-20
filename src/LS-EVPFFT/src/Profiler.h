#pragma once // Include guard to prevent multiple inclusion

#include <iostream>
#include <chrono>
#include <string>
#include <unordered_map>
#include "mpi.h"

class Profiler {
public:
  Profiler(const std::string& functionName);

  ~Profiler();

  static void print(const MPI_Comm mpiComm);

private:
  std::string functionName;
  std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
  static std::unordered_map<std::string, double> profileMap;
};

