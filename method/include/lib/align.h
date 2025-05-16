#pragma once

#include <string>
#include <vector>

#include "chrData.h"
#include "file_classes.h"

void    alignSingleWithRef(const std::string &filename, const Reference &ref, const unsigned int m,
                           const unsigned int skip_header, const unsigned int chunk_size,
                           FileMap &fileMap);
FileMap alignWithRef(const std::vector<std::string> &filenames, const Reference &ref,
                     const unsigned int m = 2, const unsigned int skip_header = 0,
                     unsigned int n_cores = 1, const unsigned int chunk_size = 64 * 1024);