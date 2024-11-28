#include <iostream>
#include <zlib.h>
#include <sstream>

#include "chrData.h"

/**
 * Parse a bed.gz file of all CpGs of a reference genome into a nested structure to be used for aligning the relevant regions of different cell files.
 *
 * @param filename path to bed.gz file to be extracted.
 * @param intervals Intervals object with search intervals.
 * @return Reference object which contain the chr information and intervals corresponding to it.
 */
Reference parseRef(const std::string &filename, Intervals intervals)
{
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file)
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
  }

  Reference ref;
  for (const auto &[chr, intervals] : intervals)
  {
    ref.emplace_back(chr, std::vector<std::vector<unsigned int>>(intervals.size()));
  }

  constexpr int bufferSize = 64 * 1024;
  char buffer[bufferSize];
  std::string partialLine;
  bool headerSkipped = false;

  unsigned int chrIndex = 0, intervalIndex = 0;
  bool done = false;

  while (!done)
  {
    int bytesRead = gzread(file, buffer, bufferSize - 1);

    if (bytesRead < 0)
    {
      std::cerr << "Error reading gzip file" << std::endl;
      gzclose(file);
    }

    buffer[bytesRead] = '\0';

    std::string fullBuffer = partialLine + buffer;
    partialLine.clear();

    std::istringstream ss(fullBuffer);
    std::string line;

    while (std::getline(ss, line))
    {
      if (ss.eof() && line.back() != '\n')
      {
        partialLine = line;
        break;
      }
      if (!headerSkipped)
      {
        headerSkipped = true;
        continue;
      }
      std::istringstream lineStream(line);
      std::string chr, temp;
      unsigned int pos;
      ss >> chr >> pos >> temp;

      if (chr != intervals[chrIndex].chr)
        continue;
      if (pos < intervals[chrIndex].intervals[intervalIndex].start)
        continue;
      while (intervals[chrIndex].intervals[intervalIndex].end <= pos)
      {
        if (chrIndex == intervals.size() - 1)
        {
          if (intervalIndex == intervals[chrIndex].intervals.size() - 1)
          {
            done = true;
            break;
          }
          intervalIndex++;
        }
        else
        {
          if (intervalIndex == intervals[chrIndex].intervals.size() - 1)
          {
            chrIndex++;
            intervalIndex = 0;
            break;
          }
          else
          {
            intervalIndex++;
          }
        }
      }

      if (done)
        break;
      if (intervals[chrIndex].intervals[intervalIndex].start <= pos && pos < intervals[chrIndex].intervals[intervalIndex].end)
      {
        ref[chrIndex].positions[intervalIndex].emplace_back(pos);
      }
    }
    if (done)
      break;
    if (bytesRead < bufferSize - 1)
      break;
  }
  return ref;
}