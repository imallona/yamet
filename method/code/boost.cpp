#include <iostream>
#include <stdexcept>

#include "boost.h"

namespace po = boost::program_options;

po::variables_map parseCommandLine(int argc, char **argv)
{
  po::variables_map vm;

  po::options_description desc("yamet - allowed options");
  desc.add_options()("help,h", "produce help message")("tsv,t", po::value<std::vector<std::string>>()->composing(), ".tsv files for different cells")("ref,r", po::value<std::string>(), "path to tsv.gz file for reference CpG sites")("bed,b", po::value<std::string>(), "path to bed file for regions of interest")("chr", po::value<std::string>(), "(optional) chr of interest")("det-out,d", po::value<std::string>(), "(optional) path to detailed output file")("out,o", po::value<std::string>(), "(optional) path to simple output file");

  po::positional_options_description p;
  p.add("tsv", -1);

  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    throw std::runtime_error("Help message displayed");
  }
  return vm;
}

std::vector<std::string> getTsvFiles(const po::variables_map &vm)
{
  return vm["tsv"].as<std::vector<std::string>>();
}

std::string getChr(const po::variables_map &vm)
{
  return vm["chr"].as<std::string>();
}

std::string getBed(const po::variables_map &vm)
{
  return vm["bed"].as<std::string>();
}

std::string getRef(const po::variables_map &vm)
{
  return vm["ref"].as<std::string>();
}

std::string getDetOut(const po::variables_map &vm)
{
  return vm["det-out"].as<std::string>();
}

std::string getOut(const po::variables_map &vm)
{
  return vm["out"].as<std::string>();
}