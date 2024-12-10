#pragma once

#include <cstdint>
#include <deque>
#include <utility>

#include "file_classes.h"

class Window {
private:
  std::deque<int8_t> container;
  unsigned int       window_size;

public:
  explicit Window(unsigned int s);

  Window &append(int8_t meth);
  void    clear();
  size_t  size();
  bool    full();
  void    notify(FileCounts &fileCounts, unsigned int binIndex, unsigned int ChrIndex);
};
