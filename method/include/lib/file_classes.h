#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <chrData.h>

struct BinCounts {
  std::vector<unsigned int> cm;
  std::vector<unsigned int> cm_1;
  unsigned long long        A        = 0;
  unsigned long long        B        = 0;
  unsigned int              m        = 0;
  unsigned int              t        = 0;
  double                    avg_meth = -1.0;
  double                    sampen   = -1.0;

  explicit BinCounts(unsigned int m) : cm(1 << m, 0), cm_1(1 << (m + 1), 0) {}
};

struct ChrCounts {
  std::string            chr;
  std::vector<BinCounts> bins;

  explicit ChrCounts(std::string c, unsigned int size, unsigned int m)
      : chr(c), bins(size, BinCounts(m)) {}
};

class FileCounts {
private:
  std::vector<ChrCounts> container;

public:
  explicit FileCounts(const Reference &ref, unsigned int m);

  void count(std::pair<unsigned int, unsigned int> idx, unsigned int chrIndex,
             unsigned int binIndex);
  void addReads(unsigned int m, unsigned int t, unsigned int chrIndex, unsigned int binIndex);
  std::vector<ChrCounts> &getContainer();
};

struct File {
  std::vector<ChrCounts> chrCounts;
  unsigned long long     A        = 0;
  unsigned long long     B        = 0;
  unsigned int           m        = 0;
  unsigned int           t        = 0;
  double                 avg_meth = -1.0;
  double                 sampen   = -1.0;

  File() {}

  explicit File(std::vector<ChrCounts> &&c) : chrCounts(std::move(c)) {}
};

using FileMapContainer = std::unordered_map<std::string, File>;

class FileMap : public FileMapContainer {
public:
  void addFile(const std::string &key, std::vector<ChrCounts> &chrCounts);
  void aggregate();
  void print(const std::vector<std::string> &filenames);
  void exportDetOut(const std::string &out, const std::vector<std::string> &filenames,
                    const Intervals &intervals);
  void exportMethOut(const std::string &out, const std::vector<std::string> &filenames,
                     const Intervals &intervals);
  void exportOut(const std::string &out, const std::vector<std::string> &filenames);
};