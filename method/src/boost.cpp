#include <boost/program_options.hpp>
#include <cerrno>
#include <cmath>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>
#include <thread>

#include "boost.h"
#include "version.h"

po::variables_map parseCommandLine(int argc, char **argv) {
  po::variables_map vm;

  po::options_description inp("input");
  // clang-format off
  inp.add_options()
    ("cell,c", po::value<std::vector<std::string>>()->composing(), 
        "tab separated files, sorted by chromosome and position, for "
        "different cells in the following format\n"
        "\n"
        " chr1    5    0    2    0\n"
        " chr1    9    1    1    1\n"
        " chr2    2    3    4    1\n"
        "\n"
        "where the columns are the chromosome, position, number of methylated reads, "
        "total number of reads and the rate respectively")
    ("reference,r", po::value<std::string>(),
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
        "where the columns are the chromosome and start position respectively")
    ("intervals,i", po::value<std::string>(),
        "bed file, sorted by chromosome and start position, for "
        "intervals of interest in the following format\n"
        "\n"
        " chr1    5     7\n"
        " chr1    10    30\n"
        " chr2    1     6\n"
        "\n"
        "where the columns are the chromosome, start position and the end position respectively")
    ("skip-header", po::value<unsigned int>()->implicit_value(1),
        "integer value indicating number of lines to skip in all file inputs"
        "(this is a default value which can be overriden by other 'skip-header-*' options)")
    ("skip-header-cell", po::value<unsigned int>()->implicit_value(1),
        "integer value indicating number of lines to skip in the cell files"
        "(overrides 'skip-header' if provided)")
    ("skip-header-reference", po::value<unsigned int>()->implicit_value(1),
        "integer value indicating number of lines to skip in the reference file"
        "(overrides 'skip-header' if provided)")
    ("skip-header-intervals", po::value<unsigned int>()->implicit_value(1),
        "integer value indicating number of lines to skip in the intervals file"
        "(overrides 'skip-header' if provided)");
  // clang-format on

  po::options_description out("output");
  // clang-format off
  out.add_options()
    ("det-out,d", po::value<std::string>(), "(optional) path to detailed output file")
    ("norm-det-out,n", po::value<std::string>(), "(optional) path to detailed normalized output file")
    ("meth-out,m", po::value<std::string>(), "(optional) path to average methylation output file")
    ("all-meth", po::value<std::string>()->default_value("false")->implicit_value("true"))
    ("out,o", po::value<std::string>(), "(optional) path to simple output file");
  // clang-format on

  po::options_description ver("verbose");
  // clang-format off
  ver.add_options()
    ("print-intervals", "print parsed intervals file")
    ("print-reference", "print parsed reference file")
    ("print-sampens", po::value<std::string>()->default_value("true")->implicit_value("true"),
        "print computed sample entropies");
  // clang-format on

  po::options_description res("resource utilisation");
  // clang-format off
  res.add_options()
    ("cores", po::value<int>()->default_value(0)->notifier(validate_cores),
        "number of cores used for simultaneously parsing methylation files")
    ("chunk-size", po::value<std::string>()->default_value("64K")->notifier(validate_chunk_size),
        "size of the buffer (per file) used for reading data. Can be specified as a "
        "positive integer (bytes) or with a suffix: B (bytes), K (kilobytes), "
        "M (megabytes), G (gigabytes). Example: 4096, 64K, 128M, 2G");
  // clang-format on

  po::options_description mis("misc");
  // clang-format off
  mis.add_options()
    ("help,h", "produce help message")
    ("version", "current version information");
  // clang-format on

  po::options_description all;
  all.add(inp).add(out).add(res).add(ver).add(mis);

  po::positional_options_description p;
  p.add("cell", -1);

  po::store(po::command_line_parser(argc, argv).options(all).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help") || argc == 1) {
    std::cout << PROJECT_NAME << std::endl << all << std::endl;
    throw std::system_error(ECANCELED, std::generic_category(), "Help message displayed");
  } else if (vm.count("version")) {
    std::cout << PROJECT_NAME << " version " << PROJECT_VERSION << std::endl
              << "created and maintained by " << PROJECT_MAINTAINER << std::endl;
    throw std::system_error(ECANCELED, std::generic_category(), "Version message displayed");
  }
  return vm;
}

// Helper get functions for different arguments

std::vector<std::string> getCellFiles(const po::variables_map &vm) {
  if (!vm.count("cell")) {
    throw std::system_error(EINVAL, std::generic_category(), "No cell coverage files provided");
  }
  return vm["cell"].as<std::vector<std::string>>();
}

std::string getRef(const po::variables_map &vm) {
  if (!vm.count("reference")) {
    throw std::system_error(EINVAL, std::generic_category(), "No reference file provided");
  }
  return vm["reference"].as<std::string>();
}

std::string getIntervals(const po::variables_map &vm) {
  if (!vm.count("intervals")) {
    throw std::system_error(EINVAL, std::generic_category(), "No intervals file provided");
  }
  return vm["intervals"].as<std::string>();
}

unsigned int getSkipHeaderTemplate(const po::variables_map &vm, const std::string &file_type) {
  if (vm.count(file_type)) {
    return vm[file_type].as<unsigned int>();
  } else if (vm.count("skip-header")) {
    return vm["skip-header"].as<unsigned int>();
  } else {
    return 0;
  }
}

unsigned int getSkipHeaderCell(const po::variables_map &vm) {
  return getSkipHeaderTemplate(vm, "skip-header-cell");
}

unsigned int getSkipHeaderReference(const po::variables_map &vm) {
  return getSkipHeaderTemplate(vm, "skip-header-reference");
}

unsigned int getSkipHeaderIntervals(const po::variables_map &vm) {
  return getSkipHeaderTemplate(vm, "skip-header-intervals");
}

std::string getDetOut(const po::variables_map &vm) {
  return vm["det-out"].as<std::string>();
}

std::string getNormDetOut(const po::variables_map &vm) {
  return vm["norm-det-out"].as<std::string>();
}

std::string getMethOut(const po::variables_map &vm) {
  return vm["meth-out"].as<std::string>();
}

bool allMeth(const po::variables_map &vm) {
  char t = (char)std::tolower(vm["all-meth"].as<std::string>()[0]);
  if (t == 't' || t == '1') {
    return true;
  } else {
    return false;
  }
}

std::string getOut(const po::variables_map &vm) {
  return vm["out"].as<std::string>();
}

unsigned int getCores(const po::variables_map &vm) {
  const unsigned int max_cores = std::thread::hardware_concurrency();
  const int          c         = vm["cores"].as<int>();
  if (c == -1) {
    return max_cores;
  } else if (c == 0) {
    return max_cores - (unsigned int)(std::floor(std::log2(max_cores)));
  } else {
    return c;
  }
}

unsigned int getChunkSize(const po::variables_map &vm) {
  static const std::regex pattern(R"(^(\d+(\.\d+)?)(B|K|M|G)?$)", std::regex::icase);
  std::smatch             match;
  const auto              chunk_size = vm["chunk-size"].as<std::string>();

  std::regex_match(chunk_size, match, pattern);

  double num  = std::stod(match[1].str());
  char   unit = match[3].str().empty() ? 'B' : std::toupper(match[3].str()[0]);

  switch (unit) {
  case 'K':
    num *= 1024;
    break;
  case 'M':
    num *= 1024 * 1024;
    break;
  case 'G':
    num *= 1024 * 1024 * 1024;
    break;
  }
  if (num < 1) {
    throw std::system_error(EINVAL, std::generic_category(),
                            "Invalid value " + chunk_size + " for 'chunk_size', less than 1 byte");
  }
  return static_cast<unsigned int>(num);
}

bool printSampens(const po::variables_map &vm) {
  char t = (char)std::tolower(vm["print-sampens"].as<std::string>()[0]);
  if (t == 't' || t == '1') {
    return true;
  } else {
    return false;
  }
}

// Validators

void validate_chunk_size(const std::string &chunk_size) {
  static const std::regex pattern(R"(^(\d+(\.\d+)?)(B|K|M|G)?$)", std::regex::icase);

  if (!std::regex_match(chunk_size, pattern)) {
    throw std::system_error(EINVAL, std::generic_category(),
                            "Invalid format " + chunk_size + " for 'chunk_size'");
  }
}

void validate_cores(const int cores) {
  const int max_cores = std::thread::hardware_concurrency();

  if (cores < -1) {
    throw std::system_error(
        EINVAL, std::generic_category(),
        "Invalid value " + std::to_string(cores) +
            " for 'cores'. Allowed values are\n  -1 : use all cores\n   0 : let the program "
            "decide\n > 0 : specify the number of cores\n");
  } else if (cores > max_cores) {
    throw std::system_error(EINVAL, std::generic_category(),
                            "Invalid value " + std::to_string(cores) +
                                " for 'cores', exceeds the maximum number of available cores, " +
                                std::to_string(max_cores));
  }
}