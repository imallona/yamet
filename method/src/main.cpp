#include <cerrno>
#include <iostream>
#include <vector>

#include "align.h"
#include "chrData.h"
#include "cli_config.h"
#include "file_classes.h"

int main(int argc, char **argv) {
  try {
    CLIConfig config = parseCommandLine(argc, argv);

    std::vector<std::string> filenames = config.cell_files;

    Intervals intervals =
        parseSearch(config.intervals, config.skip_header_intervals, config.chunk_size);
    if (config.print_intervals) {
      intervals.print(std::cout);
    }

    Reference ref =
        parseRef(config.reference, intervals, config.skip_header_reference, config.chunk_size);
    if (config.print_reference) {
      ref.print(std::cout);
    }

    ParsedInfo parsedInfo = alignWithRef(filenames, ref, 2, config.all_meth,
                                         config.skip_header_cell, config.cores, config.chunk_size);

    if (config.print_sampens || !config.det_out.empty() || !config.meth_out.empty() ||
        !config.out.empty()) {
      parsedInfo.aggregate();
      if (config.print_sampens) {
        parsedInfo.print(filenames, std::cout);
      }
      if (!config.det_out.empty()) {
        parsedInfo.exportDetOut(config.det_out, filenames, intervals);
      }
      if (!config.norm_det_out.empty()) {
        parsedInfo.exportNormDetOut(config.norm_det_out, filenames, intervals);
      }
      if (!config.meth_out.empty()) {
        parsedInfo.exportMethOut(config.meth_out, filenames, intervals);
      }
      if (!config.out.empty()) {
        parsedInfo.exportOut(config.out, filenames);
      }
    }

    return 0;
  } catch (const std::system_error &e) {
    if (e.code().value() == ECANCELED) {
      return 0;
    }
    std::cerr << "\033[31mError:\033[0m " << e.what() << std::endl;
    return 1;
  } catch (const std::runtime_error &e) {
    std::cerr << "\033[31mError:\033[0m " << e.what() << std::endl;
    return 1;
  }
}
