#include <cerrno>
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
 * @param chunk_size size of file chunk.
 * @return Reference object which contain the chr information and intervals
 * corresponding to it.
 */
Reference parseRef(const std::string &filename, const Intervals &intervals,
                   const unsigned int skip_header, const unsigned int chunk_size) {
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file) {
    throw std::system_error(errno, std::generic_category(), "Opening " + filename);
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
   * buffer size of 64K by default - this is the total amount of information from a file stored at
   * any time as a chunk
   */
  std::vector<char> buffer(chunk_size);
  std::string       fullBuffer;

  /// current chromosome being parsed from reference file
  std::string currentChr = "";
  /// indicates whether there is a region requested for a particular chromosome
  bool         firstFound      = false;
  unsigned int skipped_headers = 0;

  int          chrIndex      = -1;
  unsigned int intervalIndex = 0;
  bool         done          = false;

  while (!done) {
    /// read a chunk of data from a file
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

      /// parsing a line from reference
      std::istringstream lineStream(line);
      std::string        chr;
      unsigned int       pos;
      if (!(lineStream >> chr >> pos)) {
        throw std::system_error(EIO, std::generic_category(),
                                "in line\n\n\t\033[33m" + line +
                                    "\033[0m\n\nparsing reference file " + filename);
      }

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
    if (bytesRead < chunk_size) {
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