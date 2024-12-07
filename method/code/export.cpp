#include <fstream>
#include <iostream>

#include "chrData.h"
#include "methData.h"
#include "samp_en.h"

/**
 * Exports a tab separated file with sample entropies per region for every cell file.
 *
 * @param out path to tab separated file where detailed sample entropy outputs are to be stored.
 * @param filenames vector of filenames of all cell files.
 * @param sampens SampEns object with detailed sample entropy data.
 * @param intervals Intervals object with search intervals.
 */
void exportDetOut(const std::string &out, const std::vector<std::string> &filenames,
                  SampEns &sampens, Intervals &intervals) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    std::cerr << "Error: Could not open file " << out << " for writing." << std::endl;
    return;
  }
  outStream << "chr\tstart\tend";
  for (const auto &filename : filenames) {
    outStream << "\t" << filename;
  }
  outStream << std::endl;

  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      /// print the region information
      outStream << intervals[i].chr << "\t" << intervals[i].intervals[j].start << "\t"
                << intervals[i].intervals[j].end;
      /// print sample entropies for all files at that region
      for (const auto &filename : filenames) {
        outStream << "\t" << sampens[filename].raw[i][j];
      }
      outStream << std::endl;
    }
  }
  outStream.close();
}

/**
 * Exports a tab separated file with sample entropies aggregated for every cell file.
 *
 * @param out path to tab separated file where sample entropy outputs are to be stored.
 * @param filenames vector of filenames of all cell files.
 * @param sampens SampEns object with detailed sample entropy data.
 */
void exportOut(const std::string &out, const std::vector<std::string> &filenames,
               SampEns &sampens) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    std::cerr << "Error: Could not open file " << out << " for writing." << std::endl;
    return;
  }
  outStream << "file\tvalue" << std::endl;

  /// print aggregated sample entropy for each file
  for (const auto &filename : filenames) {
    outStream << filename << "\t" << sampens[filename].agg;
    outStream << std::endl;
  }
  outStream.close();
}