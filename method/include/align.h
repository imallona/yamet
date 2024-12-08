#pragma once

#include "chrData.h"
#include "methData.h"

void    alignSingleWithRef(const std::string &filename, Reference &ref, FileMap &fileMap);
FileMap alignWithRef(const std::vector<std::string> &filenames, Reference &ref);