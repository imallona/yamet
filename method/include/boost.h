#pragma once

#include <boost/program_options.hpp>

boost::program_options::variables_map parseCommandLine(int argc, char **argv);
std::vector<std::string>              getTsvFiles(const boost::program_options::variables_map &vm);
std::string                           getBed(const boost::program_options::variables_map &vm);
std::string                           getRef(const boost::program_options::variables_map &vm);
std::string                           getDetOut(const boost::program_options::variables_map &vm);
std::string                           getOut(const boost::program_options::variables_map &vm);