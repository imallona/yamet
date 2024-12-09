#pragma once

#include <string>
#include <vector>

#include "chrData.h"
#include "samp_en.h"

void exportDetOut(const std::string &out, const std::vector<std::string> &filenames,
                  SampEns &sampens, Intervals &intervals);

void exportOut(const std::string &out, const std::vector<std::string> &filenames, SampEns &sampens);