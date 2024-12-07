#include <future>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <zlib.h>

#include <chrData.h>
#include <methData.h>

void alignSingleWithRef(const std::string &filename, Reference &ref, FileMap &fileMap) {
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file) {
    std::cerr << "Failed to open file: " << filename << std::endl;
  }

  constexpr int bufferSize = 64 * 1024;
  char          buffer[bufferSize];
  std::string   partialLine;
  bool          headerSkipped = false;

  FileMeths meths;

  for (const auto &[chr, positions] : ref) {
    meths.emplace_back(chr, std::vector<std::vector<char>>(positions.size(), std::vector<char>()));
  }

  unsigned int chrIndex = 0, binIndex = 0, posIndex = 0;
  bool         start      = true;
  std::string  currentChr = "";
  bool         foundChr   = false;
  int          lastBinPos = -1;

  while (true) {
    int bytesRead = gzread(file, buffer, bufferSize - 1);

    if (bytesRead < 0) {
      std::cerr << "Error reading gzip file" << std::endl;
      gzclose(file);
    }

    buffer[bytesRead] = '\0';

    std::string fullBuffer = partialLine + buffer;
    partialLine.clear();

    std::istringstream ss(fullBuffer);
    std::string        line;

    while (std::getline(ss, line)) {
      if (ss.eof() && line.back() != '\n') {
        partialLine = line;
        break;
      }
      if (!headerSkipped) {
        headerSkipped = true;
        continue;
      }
      std::istringstream lineStream(line);
      std::string        chr, temp;
      unsigned int       pos;
      char               methValue;
      lineStream >> chr >> pos >> temp >> temp >> methValue;

      if (chr != currentChr) {
        binIndex = 0, posIndex = 0;
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
        // std::cout << chrIndex << "  " << binIndex << "  " << posIndex << std::endl;
        bool increment = false, exit = false;
        if (ref[chrIndex].positions[binIndex].size() > 0) {
          if (pos == ref[chrIndex].positions[binIndex][posIndex]) {
            if (lastBinPos > -1 && posIndex - lastBinPos > 1) {
              meths[chrIndex].meth[binIndex].push_back(-1);
            }
            meths[chrIndex].meth[binIndex].push_back(methValue);
            increment = true, exit = true;
            lastBinPos = posIndex;
          } else if (pos < ref[chrIndex].positions[binIndex][posIndex])
            exit = true; // break;
          else {
            increment = true;
          }
        } else {
          if (binIndex == ref[chrIndex].positions.size() - 1) {
            exit = true; // break;
          } else {
            binIndex += 1, posIndex = 0, lastBinPos = -1;
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
    if (bytesRead < bufferSize - 1 && partialLine.empty()) {
      break;
    }
  }
  gzclose(file);
  fileMap[filename] = std::move(meths);
}

FileMap alignWithRef(const std::vector<std::string> &filenames, Reference &ref) {
  FileMap fileMap;
  fileMap.reserve(filenames.size());
  std::vector<std::future<void>> futures;

  for (const auto &filename : filenames) {
    // Launch asynchronous tasks for each filename
    futures.push_back(
        std::async(std::launch::async, [&]() { alignSingleWithRef(filename, ref, fileMap); }));
  }
  for (auto &future : futures) {
    future.get();
  }
  return fileMap;
}