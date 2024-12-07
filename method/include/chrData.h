#pragma once

#include <string>
#include <vector>

struct Position {
  unsigned int start;
  unsigned int end;

  Position(unsigned int s, unsigned int e) : start(s), end(e) {}
};

struct ChrIntervals {
  std::string           chr;
  std::vector<Position> intervals;

  ChrIntervals(const std::string &c, const std::vector<Position> &i) : chr(c), intervals(i) {}
};

using Intervals = std::vector<ChrIntervals>;

struct ChrPositions {
  std::string                            chr;
  std::vector<std::vector<unsigned int>> positions;

  ChrPositions(const std::string &c, const std::vector<std::vector<unsigned int>> &p)
      : chr(c), positions(p) {}
};

using Reference = std::vector<ChrPositions>;