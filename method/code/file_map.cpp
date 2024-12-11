#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "file_classes.h"

void FileMap::addFile(const std::string &key, std::vector<ChrCounts> &chrCounts) {
  (*this)[key] = File(std::move(chrCounts));
}

void FileMap::aggregate() {
  for (auto &[filename, file] : *this) {
    for (unsigned int i = 0; i < file.chrCounts.size(); i++) {
      for (unsigned int j = 0; j < file.chrCounts[i].bins.size(); j++) {
        for (unsigned int k = 0; k < file.chrCounts[i].bins[j].cm_1.size(); k++) {
          if (k < file.chrCounts[i].bins[j].cm.size() && file.chrCounts[i].bins[j].cm[k] > 1) {
            file.chrCounts[i].bins[j].A +=
                ((unsigned long long)(file.chrCounts[i].bins[j].cm[k]) *
                 (unsigned long long)(file.chrCounts[i].bins[j].cm[k] - 1)) /
                2;
          }
          if (file.chrCounts[i].bins[j].cm_1[k] > 1) {
            file.chrCounts[i].bins[j].B +=
                ((unsigned long long)(file.chrCounts[i].bins[j].cm_1[k]) *
                 (unsigned long long)(file.chrCounts[i].bins[j].cm_1[k] - 1)) /
                2;
          }
        }
        if (file.chrCounts[i].bins[j].A > 0 && file.chrCounts[i].bins[j].B > 0) {
          file.A += file.chrCounts[i].bins[j].A;
          file.B += file.chrCounts[i].bins[j].B;
          file.chrCounts[i].bins[j].sampen =
              log(((double)file.chrCounts[i].bins[j].A) / ((double)file.chrCounts[i].bins[j].B));
        }
      }
    }
    if (file.A > 0 && file.B > 0) {
      file.sampen = log(((double)file.A) / ((double)file.B));
    }
  }
}

void FileMap::print(std::vector<std::string> filenames) {
  std::cout << "--Sample Entropies------------------" << std::endl << std::endl;
  for (const auto filename : filenames) {
    std::cout << "Filename: " << filename << std::endl;
    std::cout << "  Aggregate: " << (*this)[filename].sampen << std::endl;
    std::cout << "  Detailed:" << std::endl;
    for (unsigned int chrIndex = 0; chrIndex < (*this)[filename].chrCounts.size(); chrIndex++) {
      std::cout << "    Chromosome: " << (*this)[filename].chrCounts[chrIndex].chr << std::endl;
      for (unsigned int binIndex = 0; binIndex < (*this)[filename].chrCounts[chrIndex].bins.size();
           binIndex++) {
        std::cout << "      Bin " << binIndex << " -> "
                  << (*this)[filename].chrCounts[chrIndex].bins[binIndex].sampen << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

/**
 * Exports a tab separated file with sample entropies per region for every cell file.
 *
 * @param out path to tab separated file where detailed sample entropy outputs are to be stored.
 * @param filenames vector of filenames of all cell files.
 * @param intervals Intervals object with search intervals.
 */
void FileMap::exportDetOut(const std::string &out, const std::vector<std::string> &filenames,
                           Intervals &intervals) {
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
        outStream << "\t" << (*this)[filename].chrCounts[i].bins[j].sampen;
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
 */
void FileMap::exportOut(const std::string &out, const std::vector<std::string> &filenames) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    std::cerr << "Error: Could not open file " << out << " for writing." << std::endl;
    return;
  }
  outStream << "file\tvalue" << std::endl;

  /// print aggregated sample entropy for each file
  for (const auto &filename : filenames) {
    outStream << filename << "\t" << (*this)[filename].sampen;
    outStream << std::endl;
  }
  outStream.close();
}