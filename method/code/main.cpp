#include <iostream>

#include "align.h"
#include "boost.h"
#include "chrData.h"
#include "file_classes.h"
#include "parse_ref.h"
#include "parse_search.h"

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

    FileMap fileMap = alignWithRef(filenames, ref, 2, getNCores(vm), getNThreadsPerCore(vm));
    if (printSampens(vm)) {
      fileMap.print(filenames);
    }

    if (vm.count("det-out")) {
      fileMap.exportDetOut(getDetOut(vm), filenames, intervals);
    }
    if (vm.count("out")) {
      fileMap.exportOut(getOut(vm), filenames);
    }
    return 0;
  } catch (const po::error &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (const std::system_error &e) {
    if (e.code().value() == 99) {
      return 0;
    }
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  } catch (const std::runtime_error &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}