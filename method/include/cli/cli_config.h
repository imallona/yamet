#pragma once

#include <string>
#include <vector>

#define CLI11_ENABLE_EXTRA_VALIDATORS 1
#include <CLI/CLI.hpp>

struct CLIConfig {
  std::vector<std::string> cell_files;
  std::string              reference;
  std::string              intervals;
  unsigned int             skip_header_cell      = 0;
  unsigned int             skip_header_reference = 0;
  unsigned int             skip_header_intervals = 0;
  std::string              det_out;
  std::string              norm_det_out;
  std::string              meth_out;
  std::string              out;
  bool                     all_meth        = false;
  bool                     print_sampens   = true;
  bool                     print_intervals = false;
  bool                     print_reference = false;
  unsigned int             cores           = 0;
  unsigned int             chunk_size      = 0;
};

CLIConfig parseCommandLine(int argc, char **argv);
