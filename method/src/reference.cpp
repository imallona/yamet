#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "chrData.h"

/**
 * Parse a tab separated file of all positions of a reference genome into a nested structure
 * to be used for aligning the relevant regions of different cell files.
 *
 * @param filename path to bed.gz file to be extracted.
 * @param intervals Intervals object with search intervals.
 * @return Reference object which contain the chr information and intervals
 * corresponding to it.
 */
Reference parseRef(const std::string &filename, const Intervals &intervals) {
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file) {
    std::cerr << "Failed to open file: " << filename << std::endl;
  }

  /**
   *  create output reference variable and initialise accordingly with chromosomes where search
   * regions are present
   */
  Reference ref;
  for (const auto &[chr, chrIntervals] : intervals) {
    ref.emplace_back(chr, chrIntervals.size());
  }

  /**
   * buffer size of 64M - this is the total amount of information from a file stored at any time as
   * a chunk
   */
  constexpr int bufferSize = 64 * 1024;
  char          buffer[bufferSize];
  std::string   partialLine;

  /// current chromosome being parsed from reference file
  std::string currentChr = "";
  /// indicates whether there is a region requested for a particular chromosome
  bool firstFound = false;
  // bool headerSkipped = false;

  int          chrIndex      = -1;
  unsigned int intervalIndex = 0;
  bool         done          = false;

  while (!done) {
    /// read a chunk of data from a file
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

    /// iterate over lines in a chunk
    while (std::getline(ss, line)) {
      /// save partial line to be processed later
      if (ss.eof() && line.back() != '\n') {
        partialLine = line;
        break;
      }
      // if (!headerSkipped) {
      //   headerSkipped = true;
      //   continue;
      // }

      /// parsing a line from reference
      std::istringstream lineStream(line);
      std::string        chr, temp;
      unsigned int       pos;
      lineStream >> chr >> pos >> temp;

      /// updates chromosome currently being evaluated and resets firstFound
      if (chr != currentChr) {
        currentChr = chr;
        firstFound = false;
      }

      /**
       * once done with all positions in a chromosome, move to next chromosome with regions of
       * interest, ignoring chromosomes in between
       */
      if (chrIndex < (int)(intervals.size() - 1)) {
        if (chr == intervals[chrIndex + 1].chr) {
          chrIndex++;
          intervalIndex = 0;
          firstFound    = true;
        } else if (!firstFound) {
          continue;
        }
      } else if (chr != intervals[chrIndex].chr) {
        done = true;
        break;
      }

      /// discard positions outside regions of interest
      if (pos < intervals[chrIndex].intervals[intervalIndex].start) {
        continue;
      }

      /**
       * for a particular postion, finds the closest region on the chromosome (if it exists)
       * enclosing or before the position
       */
      while (intervals[chrIndex].intervals[intervalIndex].end <= pos) {
        if (intervalIndex < intervals[chrIndex].intervals.size() - 1) {
          intervalIndex++;
        } else {
          break;
        }
      }

      /// add valid position to output reference variable
      if (chr == intervals[chrIndex].chr &&
          intervals[chrIndex].intervals[intervalIndex].start <= pos &&
          pos < intervals[chrIndex].intervals[intervalIndex].end) {
        ref[chrIndex].positions[intervalIndex].emplace_back(pos);
      }
    }
    if (bytesRead < bufferSize - 1) {
      break;
    }
  }
  gzclose(file);
  return ref;
}

void Reference::print() {
  std::cout << "--Reference CPGs------------------" << std::endl << std::endl;
  for (const auto &[chr, positions] : *this) {
    std::cout << "Chromosome: " << chr << std::endl;
    for (unsigned int binIndex = 0; binIndex < positions.size(); binIndex++) {
      std::cout << "  Bin: " << binIndex << std::endl;
      for (const auto &pos : positions[binIndex]) {
        std::cout << "    Pos: " << pos << std::endl;
      }
    }
  }
  std::cout << std::endl;
}