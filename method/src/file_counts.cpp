#include <algorithm>
#include <vector>

#include "chrData.h"
#include "file_classes.h"

FileCounts::FileCounts(const Reference &ref, unsigned int m) {
  std::for_each(ref.begin(), ref.end(), [this, m](const auto &chrPositions) {
    container.emplace_back(chrPositions.chr, chrPositions.positions.size(), m);
  });
}

void FileCounts::count(std::pair<unsigned int, unsigned int> idx, unsigned int chrIndex,
                       unsigned int binIndex) {
  container[chrIndex].bins[binIndex].cm[idx.first] += 1;
  container[chrIndex].bins[binIndex].cm_1[idx.second] += 1;
}

void FileCounts::addReads(unsigned int m, unsigned int t, unsigned int chrIndex,
                          unsigned int binIndex) {
  container[chrIndex].bins[binIndex].m += m;
  container[chrIndex].bins[binIndex].t += t;
}

std::vector<ChrCounts> &FileCounts::getContainer() {
  return container;
}
