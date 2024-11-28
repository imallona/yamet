#include <iostream>
#include <zlib.h>
#include <sstream>

#include "datarow.h"
#include "filter.h"

void filter(const std::string &filename, std::vector<DataRow> &data, const std::string &chr)
{
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file)
  {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  constexpr int bufferSize = 32 * 1024;
  char buffer[bufferSize];
  std::string partialLine;
  bool headerSkipped = false;
  bool foundChr = false;

  while (true)
  {
    int bytesRead = gzread(file, buffer, bufferSize - 1);

    if (bytesRead < 0)
    {
      std::cerr << "Error reading gzip file" << std::endl;
      gzclose(file);
      return;
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
      DataRow row;
      std::string temp;
      lineStream >> row.chr >> row.pos >> temp >> temp >> row.rate;
      if (row.chr == chr)
      {
        foundChr = true;
        data.push_back(row);
      }
      else if (foundChr)
        break;
      else
        continue;
    }
    if (bytesRead < bufferSize - 1)
      break;
  }
  gzclose(file);
}