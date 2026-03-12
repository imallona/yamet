#include <cerrno>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "file_classes.h"
#include "file_stream.h"

/**
 * Constructs a ParsedInfo object for methylation data aggregation.
 * Initializes `clusters`, creates per-chromosome aggregation containers, and
 * reserves memory for file-level parsed results.
 *
 * @param ref Reference object
 * @param m Maximum entropy parameter or window size for aggregation calculations.
 * @param filesMeta metadata entries for all input files.
 */
ParsedInfo::ParsedInfo(const Reference &ref, unsigned int m, const FilesMeta &filesMeta)
    : filesMeta(filesMeta) {
  parseClusters(filesMeta);
  for (const auto &[chr, positions] : ref) {
    agg.emplace_back(chr, positions.size(), m, clusters.size());
  }
  fileMap.reserve(filesMeta.size());
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
  std::unordered_map<std::string, size_t> cluster_idx_map;
  cluster_idx_map.reserve(clusters.size());
  for (size_t idx = 0; idx < clusters.size(); idx++) {
    cluster_idx_map.emplace(clusters[idx], idx);
  }
  std::unordered_map<std::string, unsigned int> total_cm;
  total_cm.reserve(clusters.size());
  for (const auto &cluster : clusters) {
    total_cm.emplace(cluster, 0);
  }

  for (unsigned int i = 0; i < fileMap.begin()->second.chrCounts.size(); i++) {
    for (unsigned int j = 0; j < fileMap.begin()->second.chrCounts[i].bins.size(); j++) {
      for (const auto &cluster : clusters) {
        total_cm[cluster] = 0;
      }
      auto &agg_bin = agg[i].bins[j];
      for (auto const &fileMeta : filesMeta) {
        auto      &file     = fileMap[fileMeta.filepath];
        auto      &file_bin = file.chrCounts[i].bins[j];
        auto const cl_idx   = cluster_idx_map[fileMeta.cluster];
        for (size_t k = 0; k < file_bin.cm_1.size(); k++) {
          if (k < file_bin.cm.size()) {
            agg_bin.cm[cl_idx][k] += file_bin.cm[k];
            total_cm[fileMeta.cluster] += file_bin.cm[k];
            if (file_bin.cm[k] > 1) {
              file_bin.A += (1ULL * file_bin.cm[k] * (file_bin.cm[k] - 1)) >> 1;
            }
          }
          agg_bin.cm_1[cl_idx][k] += file_bin.cm_1[k];
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
          agg_bin.m[cl_idx] += file_bin.m;
          agg_bin.t[cl_idx] += file_bin.t;
          file_bin.avg_meth = ((double)file_bin.m) / ((double)file_bin.t);
        }
        /// expected sample entropy under independence:
        /// sampen = log(A/B) estimates -log(Pmatch), so the expectation is -log(Pmatch)
        /// where Pmatch = p_bin^2 + (1-p_bin)^2 and p_bin is the binary methylation
        /// fraction derived directly from the template counts (cm array). Using p_bin
        /// instead of avg_meth ensures that sampen_exp is calibrated to the actual
        /// binary patterns used in the sampen computation, giving a flat adjS across
        /// all average methylation levels regardless of data format or sparsity.
        {
          double total_cm_d = 0, ones_cm = 0;
          for (size_t k = 0; k < file_bin.cm.size(); k++) {
            total_cm_d += file_bin.cm[k];
            ones_cm    += __builtin_popcountll(static_cast<unsigned long long>(k)) * file_bin.cm[k];
          }
          if (total_cm_d > 0) {
            /// m = log2(cm.size()); cm.size() == 2^m, use count-trailing-zeros
            const unsigned int m_len = static_cast<unsigned int>(
                __builtin_ctzll(static_cast<unsigned long long>(file_bin.cm.size())));
            const double p_bin = ones_cm / (m_len * total_cm_d);
            const double Pm    = p_bin * p_bin + (1.0 - p_bin) * (1.0 - p_bin);
            file_bin.sampen_exp = (Pm > 0.0 && Pm < 1.0) ? -log(Pm) : 0.0;
          }
        }

        /// normalized sample entropy (additive correction):
        /// adjS = sampen - sampen_exp = log(A/B) - (-log(Pm_bin)) = log(A/B * Pm_bin)
        /// Under independence E[adjS] = 0 (flat across avg_meth regardless of sparsity)
        if (file_bin.sampen != -1 && file_bin.sampen_exp != -1) {
          file_bin.sampen_norm = file_bin.sampen - file_bin.sampen_exp;
        }

      }

      /// shannon entropy per bin across files in a cluster
      for (const auto &cluster : clusters) {
        auto const cl_idx = cluster_idx_map[cluster];
        if (total_cm[cluster] > 0) {
          agg_bin.shannon[cl_idx] = 0;
          for (size_t k = 0; k < agg_bin.cm[cl_idx].size(); k++) {
            if (((double)agg_bin.cm[cl_idx][k]) / ((double)total_cm[cluster]) > 0) {
              agg_bin.shannon[cl_idx] -=
                  (((double)agg_bin.cm[cl_idx][k]) / ((double)total_cm[cluster])) *
                  log(((double)agg_bin.cm[cl_idx][k]) / ((double)total_cm[cluster]));
            }
          }
          agg_bin.shannon[cl_idx] /= log(agg_bin.cm[cl_idx].size());
        }
        /// average methylation per bin across files
        if (agg_bin.t[cl_idx] > 0) {
          agg_bin.avg_meth[cl_idx] = ((double)agg_bin.m[cl_idx]) / ((double)agg_bin.t[cl_idx]);
        }
        /// expected shannon entropy per bin across files based on average methylation
        if (agg_bin.avg_meth[cl_idx] != -1) {
          auto &p = agg_bin.avg_meth[cl_idx];
          agg_bin.shannon_exp[cl_idx] =
              (p == 1 || p == 0) ? 0 : -2 * (p * log(p) + (1 - p) * log(1 - p));
          agg_bin.shannon_exp[cl_idx] /= log(agg_bin.cm[cl_idx].size());
        }
        /// normalized shannon entropy per bin across files by the expected shannon entropy
        if (agg_bin.shannon[cl_idx] != -1 && agg_bin.shannon_exp[cl_idx] != -1) {
          agg_bin.shannon_norm[cl_idx] =
              (agg_bin.shannon_exp[cl_idx] == 0 && agg_bin.shannon[cl_idx] == 0)
                  ? 0
                  : agg_bin.shannon[cl_idx] / agg_bin.shannon_exp[cl_idx];
        }
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

    /// expected sample entropy under independence:
    /// sampen_exp = -log(Pmatch) where Pmatch = p_bin^2 + (1-p_bin)^2, and p_bin is
    /// the binary methylation fraction derived from the aggregated template counts.
    {
      double total_cm_d = 0, ones_cm = 0;
      unsigned int m_len = 0;
      for (const auto &chrCount : file.chrCounts) {
        for (const auto &bin : chrCount.bins) {
          for (size_t k = 0; k < bin.cm.size(); k++) {
            total_cm_d += bin.cm[k];
            ones_cm    += __builtin_popcountll(static_cast<unsigned long long>(k)) * bin.cm[k];
          }
          if (m_len == 0 && !bin.cm.empty()) {
            m_len = static_cast<unsigned int>(
                __builtin_ctzll(static_cast<unsigned long long>(bin.cm.size())));
          }
        }
      }
      if (total_cm_d > 0 && m_len > 0) {
        const double p_bin = ones_cm / (m_len * total_cm_d);
        const double Pm    = p_bin * p_bin + (1.0 - p_bin) * (1.0 - p_bin);
        file.sampen_exp = (Pm > 0.0 && Pm < 1.0) ? -log(Pm) : 0.0;
      }
    }

    /// normalized sample entropy (additive correction):
    /// adjS = sampen - sampen_exp = log(A/B * Pm_bin); flat under independence
    if (file.sampen != -1 && file.sampen_exp != -1) {
      file.sampen_norm = file.sampen - file.sampen_exp;  // additive correction
    }
  }
}

/**
 * Writes a human-readable sample entropy report for each file to the provided stream.
 * This avoids direct console output so callers can capture or redirect the report as needed.
 * Outputs detailed entropy metrics including aggregate values per file and per-bin values organized
 * by chromosome.
 *
 * @param os Output stream that receives the formatted report.
 *
 * @note Output format:
 * - File-level aggregate entropy
 * - Per-chromosome breakdown
 * - Per-bin entropy values within each chromosome
 *
 * @warning Assumes all file paths in filesMeta exist in the internal fileMap.
 */
void ParsedInfo::print(std::ostream &os) {
  os << "--Sample Entropies------------------" << std::endl << std::endl;
  for (const auto &fileMeta : filesMeta) {
    os << "Filename: " << fileMeta.filepath << std::endl;
    os << "  Aggregate: " << fileMap[fileMeta.filepath].sampen << std::endl;
    os << "  Detailed:" << std::endl;
    for (unsigned int chrIndex = 0; chrIndex < fileMap[fileMeta.filepath].chrCounts.size();
         chrIndex++) {
      os << "    Chromosome: " << fileMap[fileMeta.filepath].chrCounts[chrIndex].chr << std::endl;
      for (unsigned int binIndex = 0;
           binIndex < fileMap[fileMeta.filepath].chrCounts[chrIndex].bins.size(); binIndex++) {
        os << "      Bin " << binIndex << " -> "
           << fileMap[fileMeta.filepath].chrCounts[chrIndex].bins[binIndex].sampen << std::endl;
      }
    }
    os << std::endl;
  }
}

/**
 * Exports a tab separated file which primarily provides sample entropies per region for every cell
 * file. It also provides shannon entropies and average methylation per region across cell files.
 *
 * @param out path to tab separated file where detailed sample entropy outputs are to be stored.
 * @param intervals Intervals object with search intervals.
 */
void ParsedInfo::exportDetOut(const std::string &out, const Intervals &intervals) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    throw std::runtime_error("Could not open file " + out + " for writing.");
  }
  outStream << "chr\tstart\tend";
  for (const auto &fileMeta : filesMeta) {
    outStream << "\t" << fileMeta.id;
  }
  if (clusters.size() == 1 && clusters[0].empty()) {
    outStream << "\tshannon\tavg_meth";
  } else {
    for (const auto &cluster : clusters) {
      outStream << "\tshannon." + cluster + "\tavg_meth." + cluster;
    }
  }
  outStream << std::endl;

  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      /// print the region information
      outStream << intervals[i].chr << "\t" << intervals[i].intervals[j].start << "\t"
                << intervals[i].intervals[j].end;
      /// print sample entropies for all files at that region
      for (const auto &fileMeta : filesMeta) {
        outStream << "\t" << fileMap[fileMeta.filepath].chrCounts[i].bins[j].sampen;
      }
      /// print shannon entropies and average methylation at that region
      for (size_t k = 0; k < clusters.size(); k++) {
        outStream << "\t" << agg[i].bins[j].shannon[k] << "\t" << agg[i].bins[j].avg_meth[k];
      }
      outStream << std::endl;
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
 * @param intervals Intervals object with search intervals.
 */
void ParsedInfo::exportNormDetOut(const std::string &out, const Intervals &intervals) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    throw std::runtime_error("Could not open file " + out + " for writing.");
  }
  outStream << "chr\tstart\tend";
  for (const auto &fileMeta : filesMeta) {
    outStream << "\t" << fileMeta.id;
  }
  if (clusters.size() == 1 && clusters[0].empty()) {
    outStream << "\tshannon_norm\tavg_meth";
  } else {
    for (const auto &cluster : clusters) {
      outStream << "\tshannon_norm." + cluster + "\tavg_meth." + cluster;
    }
  }
  outStream << std::endl;

  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      /// print the region information
      outStream << intervals[i].chr << "\t" << intervals[i].intervals[j].start << "\t"
                << intervals[i].intervals[j].end;
      /// print normalized sample entropies for all files at that region
      for (const auto &fileMeta : filesMeta) {
        outStream << "\t" << fileMap[fileMeta.filepath].chrCounts[i].bins[j].sampen_norm;
      }
      /// print normalized shannon entropies and average methylation at that region
      for (size_t k = 0; k < clusters.size(); k++) {
        outStream << "\t" << agg[i].bins[j].shannon_norm[k] << "\t" << agg[i].bins[j].avg_meth[k];
      }
      outStream << std::endl;
    }
  }
  outStream.close();
}

/**
 * Exports a tab separated file with average methylation per region for every cell file.
 *
 * @param out path to tab separated file where detailed sample entropy outputs are to be stored.
 * @param intervals Intervals object with search intervals.
 */
void ParsedInfo::exportMethOut(const std::string &out, const Intervals &intervals) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    throw std::runtime_error("Could not open file " + out + " for writing.");
  }
  outStream << "chr\tstart\tend";
  for (const auto &fileMeta : filesMeta) {
    outStream << "\t" << fileMeta.id;
  }
  outStream << std::endl;

  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      /// print the region information
      outStream << intervals[i].chr << "\t" << intervals[i].intervals[j].start << "\t"
                << intervals[i].intervals[j].end;
      /// print average methylation for all files at that region
      for (const auto &fileMeta : filesMeta) {
        outStream << "\t" << fileMap[fileMeta.filepath].chrCounts[i].bins[j].avg_meth;
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
 */
void ParsedInfo::exportOut(const std::string &out) {
  std::ofstream outStream(out);

  if (!outStream.is_open()) {
    throw std::runtime_error("Could not open file " + out + " for writing.");
  }
  outStream << "file\tsampen\tsampen_norm\tavg_meth" << std::endl;

  /// print aggregated sample entropy for each file
  for (const auto &fileMeta : filesMeta) {
    outStream << fileMeta.id << "\t" << fileMap[fileMeta.filepath].sampen << "\t"
              << fileMap[fileMeta.filepath].sampen_norm << "\t"
              << fileMap[fileMeta.filepath].avg_meth << std::endl;
  }
  outStream.close();
}

/**
 * Extracts unique cluster labels from metadata while preserving insertion order.
 *
 * @param filesMeta metadata entries containing cluster labels.
 */
void ParsedInfo::parseClusters(const FilesMeta &filesMeta) {
  std::unordered_set<std::string> unique_clusters;
  for (const auto &fileMeta : filesMeta) {
    if (unique_clusters.emplace(fileMeta.cluster).second) {
      clusters.emplace_back(fileMeta.cluster);
    }
  }
}
