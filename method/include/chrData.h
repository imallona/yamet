#pragma once

#include <string>
#include <vector>

struct Interval {
  unsigned int start;
  unsigned int end;

  Interval(unsigned int s, unsigned int e) : start(s), end(e) {}
};

struct ChrIntervals {
  std::string           chr;
  std::vector<Interval> intervals;

  ChrIntervals(const std::string &c, const std::vector<Interval> &i) : chr(c), intervals(i) {}
};

using Intervals = std::vector<ChrIntervals>;

struct ChrPositions {
  std::string                            chr;
  std::vector<std::vector<unsigned int>> positions;

  ChrPositions(const std::string &c, const std::vector<std::vector<unsigned int>> &p)
      : chr(c), positions(p) {}

  ChrPositions(const std::string &c, size_t size) : chr(c), positions(size) {}
};

using Reference = std::vector<ChrPositions>;