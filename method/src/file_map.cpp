#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "file_classes.h"

void ParsedInfo::addFile(const std::string &key, std::vector<ChrCounts> &chrCounts) {
  fileMap[key] = File(std::move(chrCounts));
}

/**
 * Aggregates methylation data across all files and calculates statistical metrics.
 * Performs comprehensive aggregation of methylation data including:
 * - k-mer count summation across files
 * - Sample entropy calculations
 * - Average methylation ratios
 * - Shannon entropy calculations for aggregated data
 *
 * @note This function modifies both individual file data and populates the aggregated data
 * structure
 * @warning Assumes fileMap is not empty and all files have consistent chromosome/bin structure
 *
 * @details
 * Processing steps:
 * 1. Iterates through all chromosomes and bins
 * 2. Aggregates pattern counts (cm, cm_1) across files
 * 3. Calculates A/B values for sample entropy using combinatorial formula: n*(n-1)/2
 * 4. Computes sample entropy as log(A/B) for each bin and file
 * 5. Calculates average methylation as m/t ratio
 * 6. Computes normalized Shannon entropy for aggregated pattern distributions
 */
void ParsedInfo::aggregate() {
  for (unsigned int i = 0; i < fileMap.begin()->second.chrCounts.size(); i++) {
    for (unsigned int j = 0; j < fileMap.begin()->second.chrCounts[i].bins.size(); j++) {
      unsigned int total_cm = 0;
      for (auto &[filename, file] : fileMap) {
        for (size_t k = 0; k < file.chrCounts[i].bins[j].cm_1.size(); k++) {
          if (k < file.chrCounts[i].bins[j].cm.size()) {
            agg[i].bins[j].cm[k] += file.chrCounts[i].bins[j].cm[k];
            total_cm += file.chrCounts[i].bins[j].cm[k];
            if (file.chrCounts[i].bins[j].cm[k] > 1) {
              file.chrCounts[i].bins[j].A +=
                  ((unsigned long long)(file.chrCounts[i].bins[j].cm[k]) *
                   (unsigned long long)(file.chrCounts[i].bins[j].cm[k] - 1)) /
                  2;
            }
          }
          agg[i].bins[j].cm_1[k] += file.chrCounts[i].bins[j].cm_1[k];
          if (file.chrCounts[i].bins[j].cm_1[k] > 1) {
            file.chrCounts[i].bins[j].B +=
                ((unsigned long long)(file.chrCounts[i].bins[j].cm_1[k]) *
                 (unsigned long long)(file.chrCounts[i].bins[j].cm_1[k] - 1)) /
                2;
          }
        }
        /// sample entropy per bin per file
        if (file.chrCounts[i].bins[j].A > 0 && file.chrCounts[i].bins[j].B > 0) {
          file.A += file.chrCounts[i].bins[j].A;
          file.B += file.chrCounts[i].bins[j].B;
          file.chrCounts[i].bins[j].sampen =
              log(((double)file.chrCounts[i].bins[j].A) / ((double)file.chrCounts[i].bins[j].B));
        }
        /// average methylation per bin per file
        if (file.chrCounts[i].bins[j].t > 0) {
          file.m += file.chrCounts[i].bins[j].m;
          file.t += file.chrCounts[i].bins[j].t;
          agg[i].bins[j].m += file.chrCounts[i].bins[j].m;
          agg[i].bins[j].t += file.chrCounts[i].bins[j].t;
          file.chrCounts[i].bins[j].avg_meth =
              ((double)file.chrCounts[i].bins[j].m) / ((double)file.chrCounts[i].bins[j].t);
        }
      }
      /// shannon entropy per bin across files
      if (total_cm > 0) {
        agg[i].bins[j].shannon = 0;
        for (size_t k = 0; k < agg[i].bins[j].cm.size(); k++) {
          if (((double)agg[i].bins[j].cm[k]) / ((double)total_cm) > 0) {
            agg[i].bins[j].shannon -= (((double)agg[i].bins[j].cm[k]) / ((double)total_cm)) *
                                      log(((double)agg[i].bins[j].cm[k]) / ((double)total_cm));
          }
        }
        agg[i].bins[j].shannon /= log(agg[i].bins[j].cm.size());
      }
      /// average methylation per bin across files
      if (agg[i].bins[j].t > 0) {
        agg[i].bins[j].avg_meth = ((double)agg[i].bins[j].m) / ((double)agg[i].bins[j].t);
      }
    }
  }
  for (auto &[_, file] : fileMap) {
    if (file.A > 0 && file.B > 0) {
      file.sampen = log(((double)file.A) / ((double)file.B));
    }
    if (file.t > 0) {
      file.avg_meth = ((double)file.m) / ((double)file.t);
    }
  }
}

/**
 * Prints sample entropy information for all processed files.
 * Outputs detailed entropy metrics including aggregate values per file
 * and per-bin values organized by chromosome.
 *
 * @param filenames Vector of filenames to print entropy data for
 *
 * @note Output format:
 * - File-level aggregate entropy
 * - Per-chromosome breakdown
 * - Per-bin entropy values within each chromosome
 *
 * @warning Assumes all filenames exist in the internal fileMap
 */
void ParsedInfo::print(const std::vector<std::string> &filenames) {
  std::cout << "--Sample Entropies------------------" << std::endl << std::endl;
  for (const auto &filename : filenames) {
    std::cout << "Filename: " << filename << std::endl;
    std::cout << "  Aggregate: " << fileMap[filename].sampen << std::endl;
    std::cout << "  Detailed:" << std::endl;
    for (unsigned int chrIndex = 0; chrIndex < fileMap[filename].chrCounts.size(); chrIndex++) {
      std::cout << "    Chromosome: " << fileMap[filename].chrCounts[chrIndex].chr << std::endl;
      for (unsigned int binIndex = 0; binIndex < fileMap[filename].chrCounts[chrIndex].bins.size();
           binIndex++) {
        std::cout << "      Bin " << binIndex << " -> "
                  << fileMap[filename].chrCounts[chrIndex].bins[binIndex].sampen << std::endl;
      }
    }
    std::cout << std::endl;
  }
}

/**
 * Exports a tab separated file which primarily provides sample entropies per region for every cell
 * file. It also provides shannon entropies and average methylation per region across cell files.
 *
 * @param out path to tab separated file where detailed sample entropy outputs are to be stored.
 * @param filenames vector of filenames of all cell files.
 * @param intervals Intervals object with search intervals.
 */
void ParsedInfo::exportDetOut(const std::string &out, const std::vector<std::string> &filenames,
                              const Intervals &intervals) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    std::cerr << "Error: Could not open file " << out << " for writing." << std::endl;
    return;
  }
  outStream << "chr\tstart\tend";
  for (const auto &filename : filenames) {
    outStream << "\t" << filename;
  }
  outStream << "\tshannon\tavg_meth" << std::endl;

  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      /// print the region information
      outStream << intervals[i].chr << "\t" << intervals[i].intervals[j].start << "\t"
                << intervals[i].intervals[j].end << "\t";
      /// print sample entropies for all files at that region
      for (const auto &filename : filenames) {
        outStream << fileMap[filename].chrCounts[i].bins[j].sampen << "\t";
      }
      /// print shannon entropies and average methylation at that region
      outStream << agg[i].bins[j].shannon << "\t" << agg[i].bins[j].avg_meth << std::endl;
    }
  }
  outStream.close();
}

/**
 * Exports a tab separated file with average methylation per region for every cell file.
 *
 * @param out path to tab separated file where detailed sample entropy outputs are to be stored.
 * @param filenames vector of filenames of all cell files.
 * @param intervals Intervals object with search intervals.
 */
void ParsedInfo::exportMethOut(const std::string &out, const std::vector<std::string> &filenames,
                               const Intervals &intervals) {
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
      /// print average methylation for all files at that region
      for (const auto &filename : filenames) {
        outStream << "\t" << fileMap[filename].chrCounts[i].bins[j].avg_meth;
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
void ParsedInfo::exportOut(const std::string &out, const std::vector<std::string> &filenames) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    std::cerr << "Error: Could not open file " << out << " for writing." << std::endl;
    return;
  }
  outStream << "file\tsampen\tavg_meth" << std::endl;

  /// print aggregated sample entropy for each file
  for (const auto &filename : filenames) {
    outStream << filename << "\t" << fileMap[filename].sampen << "\t" << fileMap[filename].avg_meth;
    outStream << std::endl;
  }
  outStream.close();
}