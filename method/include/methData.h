#pragma once

#include <unordered_map>
#include <vector>

struct ChrMeth {
  std::string                      chr;
  std::vector<std::vector<int8_t>> meth;

  ChrMeth(const std::string &c, const std::vector<std::vector<int8_t>> &p) : chr(c), meth(p) {}

  ChrMeth(const std::string &c, size_t size) : chr(c), meth(size) {}
};

using FileMeths = std::vector<ChrMeth>;

using FileMap = std::unordered_map<std::string, FileMeths>;