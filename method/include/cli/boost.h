#pragma once

#include <boost/program_options.hpp>
#include <string>
#include <vector>

namespace po = boost::program_options;

po::variables_map        parseCommandLine(int argc, char **argv);
std::vector<std::string> getCellFiles(const po::variables_map &vm);
std::string              getIntervals(const po::variables_map &vm);
std::string              getRef(const po::variables_map &vm);
std::string              getDetOut(const po::variables_map &vm);
std::string              getOut(const po::variables_map &vm);
unsigned int             getCores(const po::variables_map &vm);
unsigned int             getThreadsPerCore(const po::variables_map &vm);
unsigned int             getChunkSize(const po::variables_map &vm);
bool                     printSampens(const po::variables_map &vm);

void validate_chunk_size(const std::string &chunk_size);
void validate_cores(const int cores);
void validate_threads_per_core(const int n_threads_per_core);