#pragma once

#include <boost/program_options.hpp>
#include <string>
#include <vector>

namespace po = boost::program_options;

po::variables_map        parseCommandLine(int argc, char **argv);
std::vector<std::string> getTsvFiles(const po::variables_map &vm);
std::string              getBed(const po::variables_map &vm);
std::string              getRef(const po::variables_map &vm);
std::string              getDetOut(const po::variables_map &vm);
std::string              getOut(const po::variables_map &vm);
unsigned int             getNCores(const po::variables_map &vm);
unsigned int             getNThreadsPerCore(const po::variables_map &vm);
bool                     printSampens(const po::variables_map &vm);

void validate_num_cores(const int cores);
void validate_num_threads_per_core(const int n_threads_per_core);