#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <vector>

#include "chrData.h"

/**
 * Parse a bed file of search intervals into a nested structure to be used for extracting the
 * relevant regions of the reference. We require that all search regions in a chromosome be disjoint
 * from one another.
 *
 * @param filename path to bed file to be extracted.
 * @return vector of structs which contain the chr information and intervals corresponding to it.
 */
Intervals parseSearch(const std::string &filename) {
  std::ifstream bedFile(filename);
  if (!bedFile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << std::endl;
  }

  Intervals intervals;

  std::string           line;
  std::string           currentChr = "";
  std::vector<Interval> currentIntervals;

  while (std::getline(bedFile, line)) {
    /// parsing a line from regions file
    std::istringstream iss(line);
    std::string        chr;
    unsigned int       start, end;
    iss >> chr >> start >> end;

    if (chr != currentChr) {
      if (!currentIntervals.empty()) {
        intervals.emplace_back(currentChr, currentIntervals);
        currentIntervals.clear();
      }
      currentChr = chr;
    }
    currentIntervals.emplace_back(start, end);
  }

  if (!currentIntervals.empty()) {
    intervals.emplace_back(currentChr, currentIntervals);
  }

  bedFile.close();
  return intervals;
}

void Intervals::print() {
  std::cout << "--Search Regions------------------" << std::endl << std::endl;
  for (const auto &[chr, chrIntervals] : *this) {
    std::cout << "Chromosome: " << chr << std::endl;
    for (const auto &[start, end] : chrIntervals) {
      std::cout << "  Start: " << start << ", End: " << end << std::endl;
    }
  }
  std::cout << std::endl;
}