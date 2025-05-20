#include <algorithm>
#include <cerrno>
#include <iostream>
#include <mutex>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#include <align.h>
#include <chrData.h>
#include <file_classes.h>
#include <window.h>

#include "thread_pool.cpp"

/**
 * Parse a tab separated file of all covered positions of a reference genome into a nested structure
 * from which the sample entropies can be calculated. The function is called in an asynchronous
 * manner for different files.
 *
 * @param filename[in] path to tab separated file to be parsed.
 * @param ref[in] Reference object with all positions of interest in the genome.
 * @param m[in] size of window.
 * @param chunk_size[in] size of file chunk.
 * @param fileMapMutex[in] mutex to protect fileMap while writing.
 * @param fileMap[out] FileMap object of key value pairs where the parsed methylation data is stored
 * in a FileMeths object against the filename as key.
 */
void alignSingleWithRef(const std::string &filename, const Reference &ref, const unsigned int m,
                        const unsigned int skip_header, const unsigned int chunk_size,
                        std::mutex &fileMapMutex, FileMap &fileMap) {
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file) {
    throw std::system_error(errno, std::generic_category(), "Opening " + filename);
  }

  /**
   * buffer size of 64K by default - this is the total amount of information from a file stored at
   * any time as a chunk
   */
  std::vector<char> buffer(chunk_size);
  std::string       fullBuffer;

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

  while (true) {
    int bytesRead = gzread(file, buffer.data(), chunk_size);

    if (bytesRead < 0) {
      gzclose(file);
      throw std::system_error(errno, std::generic_category(), filename);
    }

    fullBuffer.append(buffer.begin(), buffer.begin() + bytesRead);

    std::istringstream ss(fullBuffer);
    std::string        line;

    /// iterate over lines in a chunk
    while (std::getline(ss, line)) {
      /// save partial line to be processed later
      if (ss.peek() == EOF) {
        if (fullBuffer.back() == '\n') {
          fullBuffer.clear();
        } else {
          fullBuffer = line;
          break;
        }
      }
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
    if (bytesRead < chunk_size) {
      break;
    }
  }
  gzclose(file);
  {
    std::lock_guard<std::mutex> lock(fileMapMutex);
    fileMap.addFile(filename, fileCounts.getContainer());
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
 * @return FileMap object of key value pairs where the parsed methylation data of every file is
 * stored in a FileMeths object against the filename as key
 */
FileMap alignWithRef(const std::vector<std::string> &filenames, const Reference &ref,
                     const unsigned int m, const unsigned int skip_header, unsigned int n_cores,
                     const unsigned int chunk_size) {
  FileMap fileMap;
  fileMap.reserve(filenames.size());

  ThreadPool                     pool(n_cores);
  std::vector<std::future<void>> results(filenames.size());
  std::mutex                     fileMapMutex;

  std::transform(
      filenames.begin(), filenames.end(), results.begin(), [&](const std::string &filename) {
        return pool.enqueue([&, filename]() {
          alignSingleWithRef(filename, ref, m, skip_header, chunk_size, fileMapMutex, fileMap);
        });
      });

  for (auto &result : results) {
    result.get();
  }

  pool.rethrow_exception(); // Check if any task threw an exception and rethrow it

  return fileMap;
}