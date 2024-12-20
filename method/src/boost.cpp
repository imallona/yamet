#include <boost/program_options.hpp>
#include <cmath>
#include <iostream>
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
    ("tsv,t", po::value<std::vector<std::string>>()->composing(), 
        "tab separated files for different cells in the following format\n"
        "\n"
        " 1    5    0    2    0\n"
        " 1    9    1    1    1\n"
        " 2    2    3    4    1\n"
        "\n"
        "where the columns are the chromosome, position, number of methylated reads, "
        "total number of reads and the rate respectively")
    ("ref,r", po::value<std::string>(),
        "tab separated file for reference sites in the following format\n"
        "\n"
        " 1    5     7\n"
        " 1    7     9\n"
        " 1    9     11\n"
        " 1    11    13\n"
        " 2    2     4\n"
        " 2    4     6\n"
        "\n"
        "where the columns are the chromosome, start position and the end position respectively")
    ("bed,b", po::value<std::string>(),
        "path to bed file for regions of interest in the following format\n"
        "\n"
        " 1    5     7\n"
        " 1    10    30\n"
        " 2    1     6\n"
        "\n"
        "where the columns are the chromosome, start position and the end position respectively");
  // clang-format on

  po::options_description out("output");
  // clang-format off
  out.add_options()
    ("det-out,d", po::value<std::string>(), "(optional) path to detailed output file")
    ("out,o", po::value<std::string>(), "(optional) path to simple output file");
  // clang-format on

  po::options_description ver("verbose");
  // clang-format off
  ver.add_options()
    ("print-bed", "print parsed regions file")
    ("print-ref", "print parsed reference file")
    ("print-sampens", po::value<std::string>()->default_value("true")->implicit_value("true"),
        "print computed sample entropies");
  // clang-format on

  po::options_description res("resource utilisation");
  // clang-format off
  res.add_options()
    ("n-cores", po::value<int>()->default_value(0)->notifier(validate_num_cores),
        "number of cores used for simultaneously parsing methylation files")
    ("n-threads-per-core",
        po::value<unsigned int>()->default_value(1)->notifier(validate_num_threads_per_core),
        "number of threads per core used for simultaneously parsing methylation files");
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
  p.add("tsv", -1);

  po::store(po::command_line_parser(argc, argv).options(all).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << PROJECT_NAME << std::endl << all << std::endl;
    std::error_code ec(99, std::generic_category());
    throw std::system_error(ec, "Help message displayed");
  } else if (vm.count("version")) {
    std::cout << PROJECT_NAME << " version " << PROJECT_VERSION << std::endl;
    std::error_code ec(199, std::generic_category());
    throw std::system_error(ec, "Version message displayed");
  } else if (!vm.count("tsv")) {
    std::error_code ec(22, std::generic_category());
    throw std::system_error(ec, "no input");
  }
  return vm;
}

std::vector<std::string> getTsvFiles(const po::variables_map &vm) {
  return vm["tsv"].as<std::vector<std::string>>();
}

std::string getBed(const po::variables_map &vm) {
  return vm["bed"].as<std::string>();
}

std::string getRef(const po::variables_map &vm) {
  return vm["ref"].as<std::string>();
}

std::string getDetOut(const po::variables_map &vm) {
  return vm["det-out"].as<std::string>();
}

std::string getOut(const po::variables_map &vm) {
  return vm["out"].as<std::string>();
}

unsigned int getNCores(const po::variables_map &vm) {
  const unsigned int max_cores = std::thread::hardware_concurrency();
  const int          c         = vm["n-cores"].as<int>();
  if (c == -1) {
    return max_cores;
  } else if (c == 0) {
    return max_cores - (unsigned int)(std::floor(std::log2(max_cores)));
  } else {
    return c;
  }
}

unsigned int getNThreadsPerCore(const po::variables_map &vm) {
  return vm["n-threads-per-core"].as<unsigned int>();
}

bool printSampens(const po::variables_map &vm) {
  char t = (char)std::tolower(vm["print-sampens"].as<std::string>()[0]);
  if (t == 't' || t == '1') {
    return true;
  } else {
    return false;
  }
}

void validate_num_cores(const int cores) {
  const int max_cores = std::thread::hardware_concurrency();

  if (cores < -1) {
    throw po::error("Error: 'n-cores' set to " + std::to_string(cores) +
                    " is invalid. It can be\n-1 : use all cores\n 0 : let the program decide\n>0 : "
                    "specified number of cores");
  } else if (cores > max_cores) {
    throw po::error("Error: 'n-cores' set to " + std::to_string(cores) +
                    " exceeds the maximum number of available cores, " + std::to_string(max_cores));
  }
}

void validate_num_threads_per_core(const int n_threads_per_core) {
  if (n_threads_per_core < 1) {
    throw po::error("Error: 'n-threads-per-core' set to " + std::to_string(n_threads_per_core) +
                    " is invalid. It must be atleast 1");
  }
}