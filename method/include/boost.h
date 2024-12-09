#pragma once

#include <boost/program_options.hpp>

boost::program_options::variables_map parseCommandLine(int argc, char **argv);
std::vector<std::string>              getTsvFiles(const boost::program_options::variables_map &vm);
std::string                           getBed(const boost::program_options::variables_map &vm);
std::string                           getRef(const boost::program_options::variables_map &vm);
std::string                           getDetOut(const boost::program_options::variables_map &vm);
std::string                           getOut(const boost::program_options::variables_map &vm);
unsigned int                          getNCores(const boost::program_options::variables_map &vm);
unsigned int getNThreadsPerCore(const boost::program_options::variables_map &vm);
bool         printSampens(const boost::program_options::variables_map &vm);

void validate_num_cores(const int cores);
void validate_num_threads_per_core(const int n_threads_per_core);