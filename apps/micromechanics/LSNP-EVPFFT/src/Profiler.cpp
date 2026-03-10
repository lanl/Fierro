#include "Profiler.h"
#include "Kokkos_Core.hpp"

// Define the static member outside the class
std::unordered_map<std::string, double> Profiler::profileMap;

Profiler::Profiler(const std::string& functionName)
  : functionName(functionName) {
#if ENABLE_PROFILING
  Kokkos::fence();
  startTime = std::chrono::high_resolution_clock::now();
#endif
}

Profiler::~Profiler() {
#if ENABLE_PROFILING
  Kokkos::fence();
  auto endTime = std::chrono::high_resolution_clock::now();
  double durationSeconds = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / 1e6; // Convert to seconds

  profileMap[functionName] += durationSeconds;
#endif
}

void Profiler::print(const MPI_Comm mpiComm) {
#if ENABLE_PROFILING
  int myRank, numRanks;
  MPI_Comm_rank(mpiComm, &myRank);
  MPI_Comm_size(mpiComm, &numRanks);

  std::map<std::string, double> sortedMap(profileMap.begin(), profileMap.end());
  std::vector<double> timeValue;
  for (const auto& pair : sortedMap) {
    timeValue.push_back(pair.second/numRanks);
  }
  MPI_Allreduce(MPI_IN_PLACE, timeValue.data(), timeValue.size(), MPI_DOUBLE, MPI_SUM, mpiComm);

  if (0 == myRank) {
    int i = 0;
    for (const auto& entry : sortedMap) {
      std::cout << entry.first << " - Time (seconds): " << timeValue[i] << std::endl;
      i++;
    }
  }
#endif
}

