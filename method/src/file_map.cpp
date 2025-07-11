#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "file_classes.h"

/**
 * Constructs a ParsedInfo object for methylation data aggregation.
 * Initializes the ParsedInfo object by initialising `agg` for each
 * chromosome in the reference and reserving memory for `fileMap`.
 *
 * @param ref Reference object
 * @param m Maximum entropy parameter or window size for aggregation calculations.
 * @param num_files Number of input files to be processed. Used for memory pre-allocation to
 * optimize performance.
 */
ParsedInfo::ParsedInfo(const Reference &ref, unsigned int m, size_t num_files) {
  for (const auto &[chr, positions] : ref) {
    agg.emplace_back(chr, positions.size(), m);
  }
  fileMap.reserve(num_files);
}

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
      auto        &agg_bin  = agg[i].bins[j];
      for (auto &[filename, file] : fileMap) {
        auto &file_bin = file.chrCounts[i].bins[j];
        for (size_t k = 0; k < file_bin.cm_1.size(); k++) {
          if (k < file_bin.cm.size()) {
            agg_bin.cm[k] += file_bin.cm[k];
            total_cm += file_bin.cm[k];
            if (file_bin.cm[k] > 1) {
              file_bin.A += (1ULL * file_bin.cm[k] * (file_bin.cm[k] - 1)) >> 1;
            }
          }
          agg_bin.cm_1[k] += file_bin.cm_1[k];
          if (file_bin.cm_1[k] > 1) {
            file_bin.B += (1ULL * file_bin.cm_1[k] * (file_bin.cm_1[k] - 1)) >> 1;
          }
        }
        /// sample entropy per bin per file
        if (file_bin.A > 0 && file_bin.B > 0) {
          file.A += file_bin.A;
          file.B += file_bin.B;
          file_bin.sampen = log(((double)file_bin.A) / ((double)file_bin.B));
        }
        /// average methylation per bin per file
        if (file_bin.t > 0) {
          file.m += file_bin.m;
          file.t += file_bin.t;
          agg_bin.m += file_bin.m;
          agg_bin.t += file_bin.t;
          file_bin.avg_meth = ((double)file_bin.m) / ((double)file_bin.t);
        }
        /// expected sample entropy per bin per file based on average methylation
        if (file_bin.avg_meth != -1) {
          auto &p             = file_bin.avg_meth;
          file_bin.sampen_exp = -log(pow(p, 2) + pow(1 - p, 2));
        }
        /// normalized sample entropy per bin per file by the expected sample entropy
        if (file_bin.sampen != -1 && file_bin.sampen_exp != -1) {
          file_bin.sampen_norm = (file_bin.sampen_exp == 0 && file_bin.sampen == 0)
                                     ? 0
                                     : file_bin.sampen / file_bin.sampen_exp;
        }
      }
      /// shannon entropy per bin across files
      if (total_cm > 0) {
        agg_bin.shannon = 0;
        for (size_t k = 0; k < agg_bin.cm.size(); k++) {
          if (((double)agg_bin.cm[k]) / ((double)total_cm) > 0) {
            agg_bin.shannon -= (((double)agg_bin.cm[k]) / ((double)total_cm)) *
                               log(((double)agg_bin.cm[k]) / ((double)total_cm));
          }
        }
        agg_bin.shannon /= log(agg_bin.cm.size());
      }
      /// average methylation per bin across files
      if (agg_bin.t > 0) {
        agg_bin.avg_meth = ((double)agg_bin.m) / ((double)agg_bin.t);
      }
      /// expected shannon entropy per bin across files based on average methylation
      if (agg_bin.avg_meth != -1) {
        auto &p             = agg_bin.avg_meth;
        agg_bin.shannon_exp = (p == 1 || p == 0) ? 0 : -2 * (p * log(p) + (1 - p) * log(1 - p));
        agg_bin.shannon_exp /= log(agg_bin.cm.size());
      }
      /// normalized shannon entropy per bin across files by the expected shannon entropy
      if (agg_bin.shannon != -1 && agg_bin.shannon_exp != -1) {
        agg_bin.shannon_norm = (agg_bin.shannon_exp == 0 && agg_bin.shannon == 0)
                                   ? 0
                                   : agg_bin.shannon / agg_bin.shannon_exp;
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
    /// expected sample entropy per file based on average methylation
    if (file.avg_meth != -1) {
      const auto &p   = file.avg_meth;
      file.sampen_exp = -log(pow(p, 2) + pow(1 - p, 2));
    }
    /// normalized sample entropy per file by the expected sample entropy
    if (file.sampen != -1 && file.sampen_exp != -1) {
      file.sampen_norm =
          (file.sampen_exp == 0 && file.sampen == 0) ? 0 : file.sampen / file.sampen_exp;
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
 * Exports a tab separated file which primarily provides normalized sample entropies per region for
 * every cell file. It also provides normalized shannon entropies and average methylation per region
 * across cell files.
 *
 * @param out path to tab separated file where detailed sample entropy outputs are to be stored.
 * @param filenames vector of filenames of all cell files.
 * @param intervals Intervals object with search intervals.
 */
void ParsedInfo::exportNormDetOut(const std::string &out, const std::vector<std::string> &filenames,
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
  outStream << "\tshannon_norm\tavg_meth" << std::endl;

  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      /// print the region information
      outStream << intervals[i].chr << "\t" << intervals[i].intervals[j].start << "\t"
                << intervals[i].intervals[j].end << "\t";
      /// print normalized sample entropies for all files at that region
      for (const auto &filename : filenames) {
        outStream << fileMap[filename].chrCounts[i].bins[j].sampen_norm << "\t";
      }
      /// print normalized shannon entropies and average methylation at that region
      outStream << agg[i].bins[j].shannon_norm << "\t" << agg[i].bins[j].avg_meth << std::endl;
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
  outStream << "file\tsampen\tsampen_norm\tavg_meth" << std::endl;

  /// print aggregated sample entropy for each file
  for (const auto &filename : filenames) {
    outStream << filename << "\t" << fileMap[filename].sampen << "\t"
              << fileMap[filename].sampen_norm << "\t" << fileMap[filename].avg_meth << std::endl;
  }
  outStream.close();
}