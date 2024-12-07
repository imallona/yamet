#include <fstream>
#include <iostream>

#include "chrData.h"
#include "methData.h"
#include "samp_en.h"

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
      outStream << intervals[i].chr << "\t" << intervals[i].intervals[j].start << "\t"
                << intervals[i].intervals[j].end;
      for (const auto &filename : filenames) {
        outStream << "\t" << sampens[filename].raw[i][j];
      }
      outStream << std::endl;
    }
  }
  outStream.close();
}

void exportOut(const std::string &out, const std::vector<std::string> &filenames,
               SampEns &sampens) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    std::cerr << "Error: Could not open file " << out << " for writing." << std::endl;
    return;
  }
  outStream << "file\tvalue" << std::endl;

  for (const auto &filename : filenames) {
    outStream << filename << "\t" << sampens[filename].agg;
    outStream << std::endl;
  }
  outStream.close();
}