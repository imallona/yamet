#pragma once

#include <vector>

#include "methData.h"
#include "datarow.h"

struct FileSampEns
{
  std::vector<std::vector<double>> raw;
  double agg;

  // FileSampEns(std::vector<std::vector<double>> s, double t) : raw(std::move(s)), agg(t) {}
};

using SampEns = std::unordered_map<std::string, FileSampEns>;

struct FileCounts
{
  unsigned int cm;
  unsigned int cm_1;
};

using Counts = std::unordered_map<std::string, FileCounts>;

SampEns sampEn(FileMap &fileMap, const unsigned int m);