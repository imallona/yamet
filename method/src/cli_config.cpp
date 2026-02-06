#include <cerrno>
#include <cmath>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <thread>

#define CLI11_ENABLE_EXTRA_VALIDATORS 1
#include <CLI/CLI.hpp>
#if (CLI11_VERSION_MAJOR >= 2 && CLI11_VERSION_MINOR >= 6)
#define HAS_CLI_READ_PERMISSIONS 1
#else
#define HAS_CLI_READ_PERMISSIONS 0
#endif

#include "cli_config.h"
#include "version.h"

auto cores_transform = CLI::Validator(
    [](std::string &input) -> std::string {
      size_t idx   = 0;
      int    cores = 0;
      try {
        cores = std::stoll(input, &idx, 10);
      } catch (const std::exception &) {
        return "Invalid value " + input + " for 'cores'. Expected an integer.";
      }
      if (idx != input.size()) {
        return "Invalid value " + input + " for 'cores'. Expected an integer.";
      }

      const unsigned int max_cores = std::thread::hardware_concurrency();
      if (cores < -1 || cores > max_cores) {
        return "Invalid value " + input +
               ". Allowed values are\n  -1 : use all cores\n   0 : let the program "
               "decide\n > 0 : specify the number of cores (up to " +
               std::to_string(max_cores) + ")";
      }

      if (cores == -1) {
        cores = max_cores;
      } else if (cores == 0) {
        cores = max_cores - static_cast<unsigned int>(std::floor(std::log2(max_cores)));
      }
      input = std::to_string(cores);
      return std::string{};
    },
    "CORES", "validate_cores");

CLIConfig parseCommandLine(int argc, char **argv) {
  CLIConfig config;
  CLI::App  app{std::string(PROJECT_NAME) + ": " + PROJECT_DESCRIPTION};

  // input options group
#if HAS_CLI_READ_PERMISSIONS
  const auto inp_validator = CLI::ExistingFile & CLI::ReadPermissions;
#else
  const auto inp_validator = CLI::ExistingFile;
#endif
  app.usage("Usage: " + std::string(PROJECT_NAME) +
            " -c <cytosine report>... -r <reference> -i <interval> [OPTIONS]");

  app.add_option("-c,--cytosine-report,--cell", config.cell_files,
                 "tab separated files, sorted by chromosome and position, for "
                 "different cells in the following format\n"
                 "\n"
                 " chr1    5    0    2    0\n"
                 " chr1    9    1    1    1\n"
                 " chr2    2    3    4    1\n"
                 "\n"
                 "where the columns are the chromosome, position, number of "
                 "methylated reads, total number of reads and the rate "
                 "respectively")
      ->required()
      ->check(inp_validator)
      ->group("Input");
  app.add_option("-r,--cytosine-locations,--reference", config.reference,
                 "tab separated file, sorted by chromosome and position, for "
                 "reference sites in the following format\n"
                 "\n"
                 " chr1    5\n"
                 " chr1    7\n"
                 " chr1    9\n"
                 " chr1    11\n"
                 " chr2    2\n"
                 " chr2    4\n"
                 "\n"
                 "where the columns are the chromosome and start position "
                 "respectively")
      ->required()
      ->check(inp_validator)
      ->group("Input");
  app.add_option("-i,--regions,--intervals", config.intervals,
                 "bed file, sorted by chromosome and start position, for "
                 "genomic regions in the following format\n"
                 "\n"
                 " chr1    5     7\n"
                 " chr1    10    30\n"
                 " chr2    1     6\n"
                 "\n"
                 "where the columns are the chromosome, start position "
                 "and the end position respectively")
      ->required()
      ->check(inp_validator)
      ->group("Input");

  unsigned int skip_header = 0;
  auto        *opt_skip_header =
      app.add_flag("--skip-header{1}", skip_header,
                   "integer value indicating number of lines to skip in all file inputs")
          ->group("Input");
  app.add_flag("--skip-header-cytosine-report{1},--skip-header-cell{1}", config.skip_header_cell,
               "integer value indicating number of lines to skip in the cell files")
      ->excludes(opt_skip_header)
      ->group("Input");
  app.add_flag("--skip-header-cytosine-locations{1},--skip-header-reference{1}",
               config.skip_header_reference,
               "integer value indicating number of lines to skip in the reference file")
      ->excludes(opt_skip_header)
      ->group("Input");
  app.add_flag("--skip-header-regions{1},--skip-header-intervals{1}", config.skip_header_intervals,
               "integer value indicating number of lines to skip in the intervals file")
      ->excludes(opt_skip_header)
      ->group("Input");

  // output options group
  app.add_option("-d,--det-out", config.det_out, "(optional) path to detailed output file")
      ->group("Output");
  app.add_option("-n,--norm-det-out", config.norm_det_out,
                 "(optional) path to detailed normalized output file")
      ->group("Output");
  app.add_option("-m,--meth-out", config.meth_out,
                 "(optional) path to average methylation output file")
      ->group("Output");
  app.add_option("-o,--out", config.out, "(optional) path to simple output file")->group("Output");
  app.add_flag("--all-meth,!--templated-meth", config.all_meth,
               "include all CpGs in methylation summaries, including those "
               "not used for template construction")
      ->group("Output");

  // verbose options group
  app.add_flag("--print-intervals", config.print_intervals, "print parsed intervals file")
      ->group("Verbose");
  app.add_flag("--print-reference", config.print_reference, "print parsed reference file")
      ->group("Verbose");
  app.add_flag("--print-sampens,!--no-print-sampens", config.print_sampens,
               "print computed sample entropies")
      ->group("Verbose");

  // resources options group
  app.add_option("--cores", config.cores,
                 "number of cores used for simultaneously parsing methylation files")
      ->transform(cores_transform)
      ->default_val(0)
      ->group("Resources");
  app.add_option("--chunk-size", config.chunk_size,
                 "size of the buffer (per file) used for reading data. Can be specified as a "
                 "positive integer (bytes) or with a suffix: B, K, M, G.\n"
                 "Example: 4096, 64K, 128M, 2G")
      ->transform(CLI::AsSizeValue(false))
      ->default_val("64K")
      ->group("Resources");

  // misc options group
  app.set_help_flag("-h,--help", "show help message")->group("Misc");
  const std::string version_info = std::string(PROJECT_NAME) + " version " + PROJECT_VERSION +
                                   "\ncreated and maintained by " + PROJECT_MAINTAINER;
  app.set_version_flag("--version", version_info, "current version information")->group("Misc");

  if (argc == 1) {
    std::cout << app.help() << std::endl;
    throw std::system_error(ECANCELED, std::generic_category(), "Help message displayed");
  }

  try {
    app.parse(argc, argv);
  } catch (const CLI::CallForHelp &) {
    std::cout << app.help() << std::endl;
    throw std::system_error(ECANCELED, std::generic_category(), "Help message displayed");
  } catch (const CLI::CallForVersion &e) {
    std::cout << e.what() << std::endl;
    throw std::system_error(ECANCELED, std::generic_category(), "Version message displayed");
  } catch (const CLI::ParseError &e) {
    throw std::runtime_error(e.what());
  }

  if (opt_skip_header->count() > 0) {
    config.skip_header_cell      = skip_header;
    config.skip_header_reference = skip_header;
    config.skip_header_intervals = skip_header;
  }
  return config;
}
