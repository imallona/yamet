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

namespace po = boost::program_options;

po::variables_map parseCommandLine(int argc, char **argv) {
  po::variables_map vm;

  po::options_description inp("input");
  // clang-format off
  inp.add_options()

    ("cytosine_report,cell,c",
      po::value<std::vector<std::string>>()->composing(),
      "Per-cell cytosine report file(s) (Bismark-like for covered cytosines). "
      "Synonyms: --cytosine_report, --cell, -c.\n"
      "Tab-separated, sorted by chromosome and position:\n"
      "  chr   pos   meth_reads   total_reads   rate")

    ("cytosine_locations,reference,r",
      po::value<std::string>(),
      "Genomic locations of all cytosines (typically CpGs). "
      "Synonyms: --cytosine_locations, --reference, -r.\n"
      "Required to reconstruct contiguous CpG sequences.\n"
      "Columns:\n"
      "  chr   pos")

    ("regions,features,target,intervals,i",
      po::value<std::string>(),
      "BED file defining genomic regions where entropies will be computed. "
      "Synonyms: --regions, --features, --target, --intervals, -i.\n"
      "Columns:\n"
      "  chr   start   end")

    ("skip-header-all,skip-header",
      po::value<unsigned int>()->implicit_value(1),
      "Header lines to skip in all input files (default: 0). "
      "Synonyms: --skip-header-all, --skip-header.")

    ("skip-header-cytosine_report,skip-header-cell",
      po::value<unsigned int>()->implicit_value(1),
      "Header lines to skip in cytosine_report/cell files (default: 0).")

    ("skip-header-cytosine_locations,skip-header-reference",
      po::value<unsigned int>()->implicit_value(1),
      "Header lines to skip in cytosine_locations/reference file (default: 0).")

    ("skip-header-regions,skip-header-features,skip-header-target,skip-header-intervals",
      po::value<unsigned int>()->implicit_value(1),
      "Header lines to skip in regions/features/target/intervals file (default: 0).");
  // clang-format on

  po::options_description out("output");
  // clang-format off
  out.add_options()
    ("det-out,d", po::value<std::string>(), "(optional) path to detailed output file")
    ("norm-det-out,n", po::value<std::string>(), "(optional) path to detailed normalized output file")
    ("meth-out,m", po::value<std::string>(), "(optional) path to average methylation output file")
    ("all-meth",
        po::value<std::string>()->default_value("false")->implicit_value("true"),
        "If true, include all CpGs in methylation summaries, including those "
        "not used for template construction (default: false).")
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
        "positive integer (bytes) or with a suffix: B, K, M, G.");
  // clang-format on

  po::options_description mis("misc");
  // clang-format off
  mis.add_options()
    ("help,h", "print help message")
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

std::vector<std::string> getCellFiles(const po::variables_map &vm) {
  if (vm.count("cell"))
    return vm["cell"].as<std::vector<std::string>>();
  if (vm.count("cytosine_report"))
    return vm["cytosine_report"].as<std::vector<std::string>>();
  throw std::system_error(EINVAL, std::generic_category(),
                          "No cell/cytosine_report files provided");
}

std::string getRef(const po::variables_map &vm) {
  if (vm.count("reference"))
    return vm["reference"].as<std::string>();
  if (vm.count("cytosine_locations"))
    return vm["cytosine_locations"].as<std::string>();
  throw std::system_error(EINVAL, std::generic_category(),
                          "No reference/cytosine_locations file provided");
}

std::string getIntervals(const po::variables_map &vm) {
  if (vm.count("intervals"))
    return vm["intervals"].as<std::string>();
  if (vm.count("regions"))
    return vm["regions"].as<std::string>();
  if (vm.count("features"))
    return vm["features"].as<std::string>();
  if (vm.count("target"))
    return vm["target"].as<std::string>();
  throw std::system_error(EINVAL, std::generic_category(),
                          "No intervals/regions/features/target file provided");
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
  return vm["skip-header-cell"].as<unsigned int>();
}

unsigned int getSkipHeaderReference(const po::variables_map &vm) {
  return vm["skip-header-reference"].as<unsigned int>();
}

unsigned int getSkipHeaderIntervals(const po::variables_map &vm) {
  return vm["skip-header-intervals"].as<unsigned int>();
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
  char t = std::tolower(vm["all-meth"].as<std::string>()[0]);
  return (t == 't' || t == '1');
}

std::string getOut(const po::variables_map &vm) {
  return vm["out"].as<std::string>();
}

unsigned int getCores(const po::variables_map &vm) {
  const unsigned int max_cores = std::thread::hardware_concurrency();
  const int c = vm["cores"].as<int>();
  if (c == -1) return max_cores;
  if (c == 0)  return max_cores - (unsigned int)(std::floor(std::log2(max_cores)));
  return c;
}

unsigned int getChunkSize(const po::variables_map &vm) {
  static const std::regex pattern(R"(^(\d+(\.\d+)?)(B|K|M|G)?$)", std::regex::icase);
  std::smatch match;
  const auto chunk_size = vm["chunk-size"].as<std::string>();

  std::regex_match(chunk_size, match, pattern);

  double num = std::stod(match[1].str());
  char unit = match[3].str().empty() ? 'B' : std::toupper(match[3].str()[0]);

  switch (unit) {
    case 'K': num *= 1024; break;
    case 'M': num *= 1024 * 1024; break;
    case 'G': num *= 1024 * 1024 * 1024; break;
  }

  if (num < 1)
    throw std::system_error(EINVAL, std::generic_category(),
                            "Invalid chunk_size < 1 byte");

  return static_cast<unsigned int>(num);
}

bool printSampens(const po::variables_map &vm) {
  char t = std::tolower(vm["print-sampens"].as<std::string>()[0]);
  return (t == 't' || t == '1');
}

void validate_chunk_size(const std::string &chunk_size) {
  static const std::regex pattern(R"(^(\d+(\.\d+)?)(B|K|M|G)?$)", std::regex::icase);
  if (!std::regex_match(chunk_size, pattern)) {
    throw std::system_error(EINVAL, std::generic_category(),
                            "Invalid format " + chunk_size + " for 'chunk_size'");
  }
}

void validate_cores(const int cores) {
  const int max_cores = std::thread::hardware_concurrency();
  if (cores < -1)
    throw std::system_error(EINVAL, std::generic_category(),
                            "Invalid cores value: " + std::to_string(cores));
  if (cores > max_cores)
    throw std::system_error(EINVAL, std::generic_category(),
                            "Requested cores exceed available hardware");
}
