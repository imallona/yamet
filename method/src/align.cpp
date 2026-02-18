#include <algorithm>
#include <cerrno>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <mutex>
#include <sstream>
#include <vector>

#include "align.h"
#include "chrData.h"
#include "file_classes.h"
#include "file_stream.h"
#include "window.h"

#include "thread_pool.cpp"

/**
 * Parse a tab separated file of all covered positions of a reference genome into a nested structure
 * from which the sample entropies can be calculated. The function is called in an asynchronous
 * manner for different files.
 *
 * @param filename[in] path to the tab separated file to be aligned with the reference and parsed
 * for k-mer counts (optionally compressed with gzip or zstd).
 * @param ref[in] Reference object with all positions of interest in the genome.
 * @param m[in] size of window.
 * @param all_meth[in] flag indicating whether positions not involved in successful templates are
 * also included.
 * @param skip_header[in] number of initial lines to skip in the input file.
 * @param chunk_size[in] size of file chunk.
 * @param fileMapMutex[in] mutex to protect fileMap while writing.
 * @param parsedInfo[out] ParsedInfo object as described in `alignWithRef`
 */
void alignSingleWithRef(const std::string &filename, const Reference &ref, const unsigned int m,
                        const bool all_meth, const unsigned int skip_header,
                        const unsigned int chunk_size, std::mutex &fileMapMutex,
                        ParsedInfo &parsedInfo) {
  FileStream file(filename, chunk_size);
  if (!file.good()) {
    throw std::system_error(errno, std::generic_category(), "Opening " + filename);
  }

  unsigned int chrIndex = 0, binIndex = 0, posIndex = 0;
  std::string  currentChr      = "";
  bool         foundChr        = false;
  unsigned int skipped_headers = 0;
  Window       window(m + 1, all_meth);

  FileCounts fileCounts(ref, m);
  /**
   * for keeping track of the last covered position discovered, used so that we can register
   * breaks between covered positions, if any
   */
  int lastBinPos = -1;

  std::string line;

  while (file.getline(line)) {
    if (skipped_headers < skip_header) {
      skipped_headers++;
      continue;
    }
    std::istringstream lineStream(line);
    std::string        chr;
    unsigned int       pos;
    unsigned int       mVal, tVal, methValue;
    if (!(lineStream >> chr >> pos >> mVal >> tVal >> methValue)) {
      throw std::system_error(EIO, std::generic_category(),
                              "in line\n\n\t\033[33m" + line + "\033[0m\n\nparsing cell file " +
                                  filename);
    }

    if (mVal > tVal) {
      throw std::system_error(
          EIO, std::generic_category(),
          "in line\n\n\t\033[33m" + line + "\033[0m\n\nof file " + filename +
              ". Total number of reads must exceed number of methylated reads.");
    }

    if (chr != currentChr) {
      binIndex = 0, posIndex = 0;
      window.clear();
      foundChr   = false;
      currentChr = chr;
    } else if (!foundChr) {
      continue;
    }

    if (!foundChr) {
      for (chrIndex = 0; chrIndex < ref.size(); chrIndex++) {
        if (ref[chrIndex].chr == chr) {
          foundChr = true;
          break;
        }
      }
      if (!foundChr) {
        continue;
      }
    }

    while (true) {
      bool increment = false, exit = false;
      if (ref[chrIndex].positions[binIndex].size() > 0) {
        if (pos == ref[chrIndex].positions[binIndex][posIndex]) {
          if (lastBinPos > -1 && posIndex - lastBinPos > 1) {
            window.append(-1, 0, 0).notify(fileCounts, chrIndex, binIndex);
          }
          window.append(methValue, mVal, tVal).notify(fileCounts, chrIndex, binIndex);
          increment = true, exit = true;
          lastBinPos = posIndex;
        } else if (pos < ref[chrIndex].positions[binIndex][posIndex])
          exit = true;
        else {
          increment = true;
        }
      } else {
        if (binIndex == ref[chrIndex].positions.size() - 1) {
          exit = true;
        } else {
          binIndex += 1, posIndex = 0, lastBinPos = -1;
          window.clear();
        }
      }
      if (increment) {
        if (binIndex == ref[chrIndex].positions.size() - 1) {
          if (posIndex == ref[chrIndex].positions[binIndex].size() - 1) {
            lastBinPos = -1, exit = true;
          } else {
            posIndex++;
          }
        } else {
          if (posIndex == ref[chrIndex].positions[binIndex].size() - 1) {
            binIndex += 1, lastBinPos = -1, posIndex = 0;
            window.clear();
          } else {
            posIndex++;
          }
        }
      }
      if (exit) {
        break;
      }
    }
  }
  file.close();
  {
    std::lock_guard<std::mutex> lock(fileMapMutex);
    parsedInfo.addFile(filename, fileCounts.getContainer());
  }
}

/**
 * A wrapper function for alignSingleWithRef which allows multiple cell files to be parsed
 * simultaneously in a multi-threaded fashion using a ThreadPool object
 *
 * @param filesMeta list of metadata entries, each containing file path, cell id and cluster
 * information.
 * @param ref Reference object with all positions of interest in the genome.
 * @param m size of window.
 * @param all_meth flag indicating whether positions not involved in successful templates are
 * also included.
 * @param skip_header number of initial lines to skip in each cell file.
 * @param n_cores maximum number of worker threads to use.
 * @param chunk_size size of file chunk.
 * @return ParsedInfo object containing two main sub-objects:
 * - fileMap, where parsed metrics for each file are stored in key-value format
 * (filename -> File object)
 * - agg, where aggregated metrics across files are stored
 * (just initialised and currently empty)
 */
ParsedInfo alignWithRef(const FilesMeta &filesMeta, const Reference &ref, const unsigned int m,
                        const bool all_meth, const unsigned int skip_header, unsigned int n_cores,
                        const unsigned int chunk_size) {
  ParsedInfo parsedInfo(ref, m, filesMeta);

  unsigned int threads = std::min(n_cores, static_cast<unsigned int>(filesMeta.size()));
  ThreadPool   pool(threads);
  std::vector<std::future<void>> results(filesMeta.size());
  std::mutex                     fileMapMutex;

  std::transform(filesMeta.begin(), filesMeta.end(), results.begin(),
                 [&](const FileMeta &fileMeta) {
                   return pool.enqueue([&, fileMeta]() {
                     alignSingleWithRef(fileMeta.filepath, ref, m, all_meth, skip_header,
                                        chunk_size, fileMapMutex, parsedInfo);
                   });
                 });

  for (auto &result : results) {
    result.get();
  }

  pool.rethrow_exception(); // Check if any task threw an exception and rethrow it

  return parsedInfo;
}

/**
 * Parses a metadata file into FilesMeta entries.
 * Each non-empty line must contain: file id, cluster id, and file path.
 * Relative file paths are resolved against the directory containing the metadata file.
 *
 * @param meta path to metadata file.
 * @param skip_header number of initial lines to skip in the metadata file.
 * @param chunk_size size of file chunk used by FileStream.
 * @return FilesMeta vector with validated and resolved file paths.
 * @throws std::system_error if metadata parsing fails or a referenced file does not exist.
 */
FilesMeta parseMeta(const std::string &meta, const unsigned int skip_header,
                    const unsigned int chunk_size) {
  FileStream file(meta, chunk_size);
  if (!file.good()) {
    throw std::system_error(errno, std::generic_category(), "Opening " + meta);
  }

  const auto meta_parent = std::filesystem::path(meta).parent_path();

  std::string line;
  FilesMeta   filesMeta{};
  size_t      line_num = 0;

  const size_t skip_header_lines = static_cast<size_t>(skip_header);

  while (file.getline(line)) {
    line_num++;
    if (line_num <= skip_header_lines) {
      continue;
    }
    /// parsing a line from metadata file
    std::istringstream iss(line);
    std::string        filepath, id, cluster;
    if (!(iss >> id >> cluster >> filepath)) {
      throw std::system_error(EIO, std::generic_category(),
                              "in line\n\n\t\033[33m" + line + "\033[0m\n\nparsing metadata file " +
                                  meta);
    }

    std::filesystem::path filepath_path(filepath);
    if (filepath_path.is_relative()) {
      filepath_path = (meta_parent / filepath_path).lexically_normal();
    }

    std::error_code ec;
    if (!std::filesystem::exists(filepath_path, ec)) {
      throw std::system_error(
          ENOENT, std::generic_category(),
          "in metadata line " + std::to_string(line_num) + "\n\n\t\033[33m" + line +
              "\033[0m\n\nreferenced file does not exist: " + filepath_path.string());
    }

    filesMeta.emplace_back(id, cluster, filepath_path.string());
  }

  file.close();
  return filesMeta;
}

/**
 * Converts a nested vector of input file paths into FilesMeta.
 * Otherwise, each outer-group index is used as a cluster identifier.
 *
 * @param files nested input where each inner vector is a cluster-specific list of files.
 * @return FilesMeta vector with id/cluster/filepath for each input file.
 */
FilesMeta parseNestVec(const std::vector<std::vector<std::string>> &files) {
  FilesMeta filesMeta{};

  if (files.size() == 1) {
    std::transform(files[0].begin(), files[0].end(), std::back_inserter(filesMeta),
                   [](const std::string &file) { return FileMeta(file); });
  } else {
    for (size_t i = 0; i < files.size(); i++) {
      std::transform(files[i].begin(), files[i].end(), std::back_inserter(filesMeta),
                     [i](const std::string &file) { return FileMeta(i, file); });
    }
  }

  return filesMeta;
}
