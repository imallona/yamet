#include <cerrno>
#include <iostream>
#include <vector>

#include "align.h"
#include "boost.h"
#include "chrData.h"
#include "file_classes.h"

int main(int argc, char **argv) {
  try {
    auto vm = parseCommandLine(argc, argv);

    std::vector<std::string> filenames = getCellFiles(vm);

    Intervals intervals =
        parseSearch(getIntervals(vm), getSkipHeaderIntervals(vm), getChunkSize(vm));
    if (vm.count("print-intervals")) {
      intervals.print();
    }

    Reference ref = parseRef(getRef(vm), intervals, getSkipHeaderReference(vm), getChunkSize(vm));
    if (vm.count("print-reference")) {
      ref.print();
    }

    ParsedInfo parsedInfo =
        alignWithRef(filenames, ref, 2, getSkipHeaderCell(vm), getCores(vm), getChunkSize(vm));

    if (printSampens(vm) || vm.count("det-out") || vm.count("meth-out") || vm.count("out")) {
      parsedInfo.aggregate();
      if (printSampens(vm)) {
        parsedInfo.print(filenames);
      }
      if (vm.count("det-out")) {
        parsedInfo.exportDetOut(getDetOut(vm), filenames, intervals);
      }
      if (vm.count("norm-det-out")) {
        parsedInfo.exportNormDetOut(getNormDetOut(vm), filenames, intervals);
      }
      if (vm.count("meth-out")) {
        parsedInfo.exportMethOut(getMethOut(vm), filenames, intervals);
      }
      if (vm.count("out")) {
        parsedInfo.exportOut(getOut(vm), filenames);
      }
    }

    return 0;
  } catch (const std::system_error &e) {
    if (e.code().value() == ECANCELED) {
      return 0;
    }
    std::cerr << "\033[31mError:\033[0m " << e.what() << std::endl;
    return 1;
  } catch (const std::runtime_error &e) {
    std::cerr << "\033[31mError:\033[0m " << e.what() << std::endl;
    return 1;
  }
}