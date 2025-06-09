#include <cerrno>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "chrData.h"

#include "file_stream.cpp"

/**
 * Parse a bed file of search intervals into a nested structure to be used for extracting the
 * relevant regions of the reference. We require that all search regions in a chromosome be disjoint
 * from one another.
 *
 * @param filename path to bed file to be extracted (optionally compressed with gzip or zstd).
 * @return vector of structs which contain the chr information and intervals corresponding to it.
 */
Intervals parseSearch(const std::string &filename, const unsigned int skip_header,
                      const unsigned int chunk_size) {
  FileStream bedFile(filename, chunk_size);
  if (!bedFile.good()) {
    throw std::system_error(errno, std::generic_category(), "Opening " + filename);
  }

  Intervals intervals;

  std::string           line;
  std::string           currentChr = "";
  int                   lastStart = -1, lastEnd = -1;
  unsigned int          skipped_headers = 0;
  std::vector<Interval> currentIntervals;

  while (bedFile.getline(line)) {
    if (skipped_headers < skip_header) {
      skipped_headers++;
      continue;
    }
    /// parsing a line from regions file
    std::istringstream iss(line);
    std::string        chr;
    unsigned int       start, end;
    if (!(iss >> chr >> start >> end)) {
      throw std::system_error(EIO, std::generic_category(),
                              "in line\n\n\t\033[33m" + line +
                                  "\033[0m\n\nparsing intervals file " + filename);
    }

    if (start >= end) {
      throw std::system_error(EIO, std::generic_category(),
                              "in line\n\n\t\033[33m" + line + "\033[0m\n\nof file " + filename +
                                  ". Interval starts after interval end.");
    }

    if (chr != currentChr) {
      if (!currentIntervals.empty()) {
        intervals.emplace_back(currentChr, currentIntervals);
        currentIntervals.clear();
      }
      currentChr = chr;
      lastStart = lastEnd = -1;
    }
    if (lastStart != -1 && lastEnd != -1 && lastEnd > start) {
      throw std::system_error(EIO, std::generic_category(),
                              "in line\n\n\t\033[33m" + line + "\033[0m\n\nof file " + filename +
                                  ". Overlaps with previous interval.");
    }
    lastStart = start, lastEnd = end;
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