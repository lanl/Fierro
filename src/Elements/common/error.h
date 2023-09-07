#pragma once

#include "common.h"

#include <cassert>
#include <exception>

class FactorizationError : public std::runtime_error {
 public:
  FactorizationError(const std::string info) : std::runtime_error(info) {}
};

class SolutionError : public std::runtime_error {
 public:
  SolutionError(const std::string info) : std::runtime_error(info) {}
};

class VtkGridError : public std::runtime_error {
 public:
  VtkGridError(const std::string info) : std::runtime_error(info) {}
};
