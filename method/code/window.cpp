#include <cstdint>
#include <deque>
#include <utility>

#include <file_classes.h>
#include <window.h>

Window::Window(unsigned int s) : window_size(s) {}

Window &Window::append(int8_t meth) {
  if (container.size() >= window_size) {
    container.pop_front();
  }
  container.push_back(meth);
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

void Window::notify(FileCounts &fileCounts, unsigned int chrIndex, unsigned int binIndex) {
  if (container.size() < window_size) {
    return;
  }

  std::pair<unsigned int, unsigned int> idx = {0, 0};
  bool                                  add = true;

  for (unsigned int k; k < window_size - 1; k++) {
    if (container[k] == -1) {
      add = false;
      break;
    }
    idx.first += (container[k]) * (1 << k);
  }
  if (add && container[window_size - 1] != -1) {
    idx.second = idx.first + (container[window_size - 1]) * (1 << (window_size - 1));
  } else {
    return;
  }

  fileCounts.count(idx, chrIndex, binIndex);
}