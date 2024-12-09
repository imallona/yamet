#pragma once

#include <string>
#include <vector>

#include "chrData.h"
#include "methData.h"

void    alignSingleWithRef(const std::string &filename, Reference &ref, FileMap &fileMap);
FileMap alignWithRef(const std::vector<std::string> &filenames, Reference &ref,
                     unsigned int n_cores, unsigned int n_threads_per_core);