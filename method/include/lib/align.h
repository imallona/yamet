#pragma once

#include <string>
#include <vector>

#include "chrData.h"
#include "file_classes.h"

void    alignSingleWithRef(const std::string &filename, const Reference &ref, const unsigned int m,
                           const unsigned int chunk_size, FileMap &fileMap);
FileMap alignWithRef(const std::vector<std::string> &filenames, const Reference &ref,
                     const unsigned int m, unsigned int n_cores, unsigned int n_threads_per_core,
                     const unsigned int chunk_size);