#pragma once

#include <cstdint>
#include <deque>
#include <utility>

#include "file_classes.h"

struct BinMeth {
  int8_t       meth;
  unsigned int m;
  unsigned int t;

  explicit BinMeth(int8_t meth, unsigned int m, unsigned int t) : meth(meth), m(m), t(t) {}
};

class Window {
private:
  std::deque<BinMeth> container;
  unsigned int        window_size;
  bool                all_meth;
  void (Window::*notifyPointer)(FileCounts &, unsigned int, unsigned int);

  bool notifyCommon(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex);
  void notifyAllMeth(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex);
  void notifyDefault(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex);

public:
  explicit Window(unsigned int s, bool all_meth);

  Window &append(int8_t meth, unsigned int m, unsigned int t);
  void    clear();
  size_t  size();
  bool    full();
  void    notify(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex);
};
