#include <iostream>

#include "align.h"
#include "boost.h"
#include "chrData.h"
#include "file_classes.h"

int main(int argc, char **argv) {
  try {
    auto vm = parseCommandLine(argc, argv);

    std::vector<std::string> filenames = getTsvFiles(vm);

    Intervals intervals = parseSearch(getBed(vm));
    if (vm.count("print-bed")) {
      intervals.print();
    }

    Reference ref = parseRef(getRef(vm), intervals);
    if (vm.count("print-ref")) {
      ref.print();
    }

    FileMap fileMap = alignWithRef(filenames, ref, 2, getNCores(vm), getNThreadsPerCore(vm));

    if (printSampens(vm) || vm.count("det-out") || vm.count("out")) {
      fileMap.aggregate();
      if (printSampens(vm)) {
        fileMap.print(filenames);
      }
      if (vm.count("det-out")) {
        fileMap.exportDetOut(getDetOut(vm), filenames, intervals);
      }
      if (vm.count("out")) {
        fileMap.exportOut(getOut(vm), filenames);
      }
    }

    return 0;
  } catch (const po::error &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  } catch (const std::system_error &e) {
    if (e.code().value() == 99 || e.code().value() == 199) {
      return 0;
    } else if (e.code().value() == 499) {
      std::cout << "Error: no input!" << std::endl;
      return 1;
    }
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  } catch (const std::runtime_error &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}