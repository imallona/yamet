#pragma once

#include <boost/program_options.hpp>
#include <string>
#include <vector>

namespace po = boost::program_options;

po::variables_map        parseCommandLine(int argc, char **argv);
std::vector<std::string> getCellFiles(const po::variables_map &vm);
std::string              getRef(const po::variables_map &vm);
std::string              getIntervals(const po::variables_map &vm);
unsigned int getSkipHeaderTemplate(const po::variables_map &vm, const std::string &file_type);
unsigned int getSkipHeaderCell(const po::variables_map &vm);
unsigned int getSkipHeaderReference(const po::variables_map &vm);
unsigned int getSkipHeaderIntervals(const po::variables_map &vm);
std::string  getDetOut(const po::variables_map &vm);
std::string  getNormDetOut(const po::variables_map &vm);
std::string  getMethOut(const po::variables_map &vm);
bool         allMeth(const po::variables_map &vm);
std::string  getOut(const po::variables_map &vm);
unsigned int getCores(const po::variables_map &vm);
unsigned int getChunkSize(const po::variables_map &vm);
bool         printSampens(const po::variables_map &vm);

void validate_chunk_size(const std::string &chunk_size);
void validate_cores(const int cores);