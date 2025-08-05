#include <cstdint>
#include <deque>
#include <iostream>
#include <utility>

#include "file_classes.h"
#include "window.h"

Window::Window(unsigned int s, bool a) : window_size(s), all_meth(a) {
  notifyPointer = all_meth ? &Window::notifyAllMeth : &Window::notifyDefault;
}

Window &Window::append(int8_t meth, unsigned int m, unsigned int t) {
  if (container.size() >= window_size) {
    container.pop_front();
  }
  container.emplace_back(meth, m, t);
  return *this;
}

void Window::clear() {
  container.clear();
}

size_t Window::size() {
  return container.size();
}

bool Window::full() {
  return container.size() == window_size;
}

bool Window::notifyCommon(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex) {
  if (!full()) {
    return false;
  }

  std::pair<unsigned int, unsigned int> idx = {0, 0};

  for (unsigned int k = 0; k < window_size - 1; k++) {
    if (container[k].meth == -1) {
      return false;
    }
    idx.first += (container[k].meth) * (1 << k);
  }
  if (container[window_size - 1].meth == -1) {
    return false;
  }
  idx.second = idx.first + (container[window_size - 1].meth) * (1 << (window_size - 1));
  fileCounts.count(idx, chrIndex, binIndex);
  return true;
}

void Window::notifyAllMeth(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex) {
  fileCounts.addReads(container.back().m, container.back().t, chrIndex, binIndex);
  notifyCommon(fileCounts, chrIndex, binIndex);
}

void Window::notifyDefault(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex) {
  if (notifyCommon(fileCounts, chrIndex, binIndex)) {
    fileCounts.addReads(container.front().m, container.front().t, chrIndex, binIndex);
  }
}

void Window::notify(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex) {
  (this->*notifyPointer)(fileCounts, chrIndex, binIndex);
}