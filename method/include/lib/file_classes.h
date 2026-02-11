#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "chrData.h"

struct BinCounts {
  std::vector<unsigned int> cm;
  std::vector<unsigned int> cm_1;
  unsigned long long        A           = 0;
  unsigned long long        B           = 0;
  unsigned int              m           = 0;
  unsigned int              t           = 0;
  double                    avg_meth    = -1.0;
  double                    sampen      = -1.0;
  double                    sampen_exp  = -1.0;
  double                    sampen_norm = -1.0;

  explicit BinCounts(unsigned int m) : cm(1 << m, 0), cm_1(1 << (m + 1), 0) {}
};

struct ChrCounts {
  std::string            chr;
  std::vector<BinCounts> bins;

  explicit ChrCounts(const std::string &c, unsigned int size, unsigned int m)
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
  unsigned long long     A           = 0;
  unsigned long long     B           = 0;
  unsigned int           m           = 0;
  unsigned int           t           = 0;
  double                 avg_meth    = -1.0;
  double                 sampen      = -1.0;
  double                 sampen_exp  = -1.0;
  double                 sampen_norm = -1.0;

  File() {}

  explicit File(std::vector<ChrCounts> &&c) : chrCounts(std::move(c)) {}
};

using FileMap = std::unordered_map<std::string, File>;

struct BinAgg {
  std::vector<std::vector<unsigned int>> cm;
  std::vector<std::vector<unsigned int>> cm_1;
  std::vector<unsigned int>              m;
  std::vector<unsigned int>              t;
  std::vector<double>                    avg_meth;
  std::vector<double>                    shannon;
  std::vector<double>                    shannon_exp;
  std::vector<double>                    shannon_norm;

  explicit BinAgg(const unsigned int m, const unsigned int n_clusters)
      : cm(n_clusters, std::vector<unsigned int>(1 << m, 0)),
        cm_1(n_clusters, std::vector<unsigned int>(1 << (m + 1), 0)), m(n_clusters, 0),
        t(n_clusters, 0), avg_meth(n_clusters, -1.0), shannon(n_clusters, -1.0),
        shannon_exp(n_clusters, -1.0), shannon_norm(n_clusters, -1.0) {}
};

struct ChrAgg {
  std::string         chr;
  std::vector<BinAgg> bins;

  explicit ChrAgg(const std::string &c, const unsigned int size, const unsigned int m,
                  const unsigned int n_clusters)
      : chr(c), bins(size, BinAgg(m, n_clusters)) {}
};

struct FileMeta {
  std::string id;
  std::string cluster;
  std::string filepath;

  explicit FileMeta(const std::string &i, const std::string &c, const std::string &p)
      : id(i), cluster(c), filepath(p) {}

  explicit FileMeta(size_t c, const std::string &p)
      : cluster(std::to_string(c)), id(p), filepath(p) {}

  explicit FileMeta(const std::string &p) : id(p), filepath(p) {}
};

using FilesMeta = std::vector<FileMeta>;

class ParsedInfo {
public:
  FileMap                  fileMap;
  const FilesMeta         &filesMeta;
  std::vector<std::string> clusters;
  std::vector<ChrAgg>      agg;

  explicit ParsedInfo(const Reference &ref, unsigned int m, const FilesMeta &filesMeta);

  void addFile(const std::string &key, std::vector<ChrCounts> &chrCounts);
  void aggregate();
  void print(std::ostream &os);
  void exportDetOut(const std::string &out, const Intervals &intervals);
  void exportNormDetOut(const std::string &out, const Intervals &intervals);
  void exportMethOut(const std::string &out, const Intervals &intervals);
  void exportOut(const std::string &out);

private:
  void parseClusters(const FilesMeta &filesMeta);
};
