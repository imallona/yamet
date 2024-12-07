#include <iostream>

#include "align.h"
#include "boost.h"
#include "chrData.h"
#include "export.h"
#include "methData.h"
#include "parse_ref.h"
#include "parse_search.h"
#include "samp_en.h"
// #include "filter.h"

int main(int argc, char **argv) {
  try {
    auto vm = parseCommandLine(argc, argv);

    std::vector<std::string> filenames = getTsvFiles(vm);

    Intervals intervals = parseSearch(getBed(vm));
    if (vm.count("print-bed")) {
      std::cout << "--Search Regions------------------" << std::endl << std::endl;
      for (const auto &[chr, intervals] : intervals) {
        std::cout << "Chromosome: " << chr << std::endl;
        for (const auto &[start, end] : intervals) {
          std::cout << "  Start: " << start << ", End: " << end << std::endl;
        }
      }
      std::cout << std::endl;
    }

    Reference ref = parseRef(getRef(vm), intervals);
    if (vm.count("print-ref")) {
      std::cout << "--Reference CPGs------------------" << std::endl << std::endl;
      for (const auto &[chr, positions] : ref) {
        std::cout << "Chromosome: " << chr << std::endl;
        for (unsigned int binIndex = 0; binIndex < positions.size(); binIndex++) {
          std::cout << "  Bin: " << binIndex << std::endl;
          for (const auto &pos : positions[binIndex]) {
            std::cout << "    Pos: " << pos << std::endl;
          }
        }
      }
      std::cout << std::endl;
    }

    FileMap fileMap = alignWithRef(filenames, ref);
    if (vm.count("print-tsv")) {
      std::cout << "--Cell Files------------------" << std::endl << std::endl;
      for (const auto &filename : filenames) {
        std::cout << "Filename: " << filename << std::endl;
        for (unsigned int chrIndex = 0; chrIndex < fileMap[filename].size(); chrIndex++) {
          std::cout << "  Chromosome: " << fileMap[filename][chrIndex].chr << std::endl;
          for (unsigned int binIndex = 0; binIndex < fileMap[filename][chrIndex].meth.size();
               binIndex++) {
            std::cout << "    Bin: " << binIndex
                      << ", size: " << fileMap[filename][chrIndex].meth[binIndex].size()
                      << std::endl;
            for (const auto &methValue : fileMap[filename][chrIndex].meth[binIndex]) {
              std::cout << "      Meth: " << methValue << std::endl;
            }
          }
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    SampEns sampens = sampEn(fileMap, 2);
    if (printSampens(vm)) {
      std::cout << "--Sample Entropies------------------" << std::endl << std::endl;
      for (const auto &file : fileMap) {
        std::cout << "Filename: " << file.first << std::endl;
        std::cout << "  Aggregate: " << sampens[file.first].agg << std::endl;
        std::cout << "  Detailed:" << std::endl;
        for (unsigned int chrIndex = 0; chrIndex < sampens[file.first].raw.size(); chrIndex++) {
          std::cout << "    Chromosome: " << file.second[chrIndex].chr << std::endl;
          for (unsigned int binIndex = 0; binIndex < sampens[file.first].raw[chrIndex].size();
               binIndex++) {
            std::cout << "      Bin " << binIndex << " -> "
                      << sampens[file.first].raw[chrIndex][binIndex] << std::endl;
          }
        }
        std::cout << std::endl;
      }
    }

    if (vm.count("det-out")) {
      exportDetOut(getDetOut(vm), filenames, sampens, intervals);
    }
    if (vm.count("out")) {
      exportOut(getOut(vm), filenames, sampens);
    }
    return 0;
  } catch (const std::runtime_error &e) {
    return 0;
  }
}