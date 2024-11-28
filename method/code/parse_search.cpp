#include <iostream>
#include <fstream>
#include <sstream>
#include <list>

#include "chrData.h"

/**
 * Parse a bed file of search intervals into a nested structure to be used for extracting the relevant regions of the reference.
 *
 * @param filename path to bed file to be extracted.
 * @return vector of structs which contain the chr information and intervals corresponding to it.
 */
std::vector<ChrIntervals> parseSearch(const std::string &filename)
{
  std::ifstream bedFile(filename);
  if (!bedFile.is_open())
  {
    std::cerr << "Error: Could not open file " << filename << std::endl;
  }

  std::vector<ChrIntervals> intervals;

  std::string line;
  std::string currentChr = "";
  std::vector<Position> currentIntervals;

  while (std::getline(bedFile, line))
  {
    std::istringstream iss(line);
    std::string chr;
    unsigned int start, end;
    iss >> chr >> start >> end;

    if (chr != currentChr)
    {
      if (!currentIntervals.empty())
      {
        intervals.emplace_back(currentChr, currentIntervals);
        currentIntervals.clear();
      }
      currentChr = chr;
    }
    currentIntervals.emplace_back(start, end);
  }

  if (!currentIntervals.empty())
  {
    intervals.emplace_back(currentChr, currentIntervals);
  }

  bedFile.close();
  return intervals;
}
