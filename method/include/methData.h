#pragma once

#include <vector>
#include <unordered_map>

struct ChrMeth
{
  std::string chr;
  std::vector<std::vector<char>> meth;

  ChrMeth(const std::string &c, const std::vector<std::vector<char>> &p) : chr(c), meth(p) {}
};

using FileMeths = std::vector<ChrMeth>;

using FileMap = std::unordered_map<std::string, FileMeths>;