#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <chrData.h>

struct BinCounts {
  std::vector<unsigned int> cm;
  std::vector<unsigned int> cm_1;
  unsigned long long        A;
  unsigned long long        B;
  unsigned int              m;
  unsigned int              t;
  double                    sampen;

  BinCounts(unsigned int m)
      : cm(1 << m, 0), cm_1(1 << (m + 1), 0), A(0), B(0), m(0), t(0), sampen(-1.0) {}
};

struct ChrCounts {
  std::string            chr;
  std::vector<BinCounts> bins;

  ChrCounts(std::string c, unsigned int size, unsigned int m) : chr(c), bins(size, BinCounts(m)) {}
};

class FileCounts {
private:
  std::vector<ChrCounts> container;

public:
  explicit FileCounts(Reference &ref, unsigned int m);

  void count(std::pair<unsigned int, unsigned int> idx, unsigned int chrIndex,
             unsigned int binIndex);
  void addReads(unsigned int m, unsigned int t, unsigned int chrIndex, unsigned int binIndex);
  std::vector<ChrCounts> &getContainer();
};

struct File {
  std::vector<ChrCounts> chrCounts;
  unsigned long long     A;
  unsigned long long     B;
  double                 sampen;

  File() : A(0), B(0), sampen(-1.0) {}

  File(std::vector<ChrCounts> &&c) : chrCounts(std::move(c)), A(0), B(0), sampen(-1.0) {}
};

using FileMapContainer = std::unordered_map<std::string, File>;

class FileMap : public FileMapContainer {
public:
  void addFile(const std::string &key, std::vector<ChrCounts> &chrCounts);
  void aggregate();
  void print(std::vector<std::string> filenames);
  void exportDetOut(const std::string &out, const std::vector<std::string> &filenames,
                    Intervals &intervals);
  void exportShannon(const std::string &out, Intervals &intervals);
  void exportOut(const std::string &out, const std::vector<std::string> &filenames);
};