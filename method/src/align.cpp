#include <algorithm>
#include <cerrno>
#include <iostream>
#include <mutex>
#include <sstream>
#include <vector>

#include <align.h>
#include <chrData.h>
#include <file_classes.h>
#include <window.h>

#include "file_stream.cpp"
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
 * @param chunk_size[in] size of file chunk.
 * @param fileMapMutex[in] mutex to protect fileMap while writing.
 * @param parsedInfo[out] ParsedInfo object as described in `alignWithRef`
 */
void alignSingleWithRef(const std::string &filename, const Reference &ref, const unsigned int m,
                        const unsigned int skip_header, const unsigned int chunk_size,
                        std::mutex &fileMapMutex, ParsedInfo &parsedInfo) {
  FileStream file(filename, chunk_size);
  if (!file.good()) {
    throw std::system_error(errno, std::generic_category(), "Opening " + filename);
  }

  unsigned int chrIndex = 0, binIndex = 0, posIndex = 0;
  std::string  currentChr      = "";
  bool         foundChr        = false;
  unsigned int skipped_headers = 0;
  Window       window(m + 1);

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
            window.append(-1).notify(fileCounts, chrIndex, binIndex);
          }
          window.append(methValue).notify(fileCounts, chrIndex, binIndex);
          fileCounts.addReads(mVal, tVal, chrIndex, binIndex);
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
 * @param filenames list of tab separated files to be parsed.
 * @param ref Reference object with all positions of interest in the genome.
 * @param m size of window.
 * @param chunk_size size of file chunk.
 * @return ParsedInfo object containing two main sub-objects:
 * - fileMap, where parsed metrics for each file are stored in key-value format
 * (filename -> File object)
 * - agg, where aggregated metrics across files are stored
 * (just initialised and currently empty)
 */
ParsedInfo alignWithRef(const std::vector<std::string> &filenames, const Reference &ref,
                        const unsigned int m, const unsigned int skip_header, unsigned int n_cores,
                        const unsigned int chunk_size) {
  ParsedInfo parsedInfo(ref, m, filenames.size());

  unsigned int threads = std::min(n_cores, static_cast<unsigned int>(filenames.size()));
  ThreadPool   pool(threads);
  std::vector<std::future<void>> results(filenames.size());
  std::mutex                     fileMapMutex;

  std::transform(
      filenames.begin(), filenames.end(), results.begin(), [&](const std::string &filename) {
        return pool.enqueue([&, filename]() {
          alignSingleWithRef(filename, ref, m, skip_header, chunk_size, fileMapMutex, parsedInfo);
        });
      });

  for (auto &result : results) {
    result.get();
  }

  pool.rethrow_exception(); // Check if any task threw an exception and rethrow it

  return parsedInfo;
}