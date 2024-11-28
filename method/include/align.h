#pragma once

#include "methData.h"
#include "chrData.h"

void alignSingleWithRef(const std::string &filename, std::vector<ChrPositions> &positions, std::unordered_map<std::string, std::vector<ChrMeth>> &fileMap);
FileMap alignWithRef(const std::vector<std::string> &filenames, Reference &ref);