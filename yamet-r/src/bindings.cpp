#include <Rcpp.h>
#include <vector>

#include "yamet/align.h"
#include "yamet/chrData.h"
#include "yamet/file_classes.h"

// [[Rcpp::export]]
Rcpp::List yamet_cpp(const Rcpp::CharacterVector filenames, const std::string &reference_path,
                     const std::string &intervals_path, const unsigned int cores,
                     const unsigned int chunk_size) {
  std::vector<std::string> cpp_filenames = Rcpp::as<std::vector<std::string>>(filenames);
  Intervals                intervals     = parseSearch(intervals_path, 0);
  Reference                ref           = parseRef(reference_path, intervals, 0, chunk_size);
  FileMap                  fileMap = alignWithRef(cpp_filenames, ref, 2, 0, cores, chunk_size);

  fileMap.aggregate();

  unsigned int nrows = 0;
  for (unsigned int i = 0; i < intervals.size(); i++) {
    nrows += intervals[i].intervals.size();
  }

  Rcpp::CharacterVector chr(nrows);
  Rcpp::IntegerVector   start(nrows);
  Rcpp::IntegerVector   end(nrows);
  Rcpp::NumericMatrix   sampens(nrows, filenames.size());
  Rcpp::NumericMatrix   meths(nrows, filenames.size());

  unsigned int row = 0;
  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      for (unsigned int k = 0; k < cpp_filenames.size(); k++) {
        sampens(row, k) = fileMap[cpp_filenames[k]].chrCounts[i].bins[j].sampen;
        if (sampens(row, k) == -1) {
          sampens(row, k) = NA_REAL;
        }
        meths(row, k) = fileMap[cpp_filenames[k]].chrCounts[i].bins[j].avg_meth;
        if (meths(row, k) == -1) {
          meths(row, k) = NA_REAL;
        }
      }
      chr[row]   = intervals[i].chr;
      start[row] = intervals[i].intervals[j].start;
      end[row]   = intervals[i].intervals[j].end;
      row++;
    }
  }

  return Rcpp::List::create(Rcpp::Named("chr") = chr, Rcpp::Named("start") = start,
                            Rcpp::Named("end") = end, Rcpp::Named("sampens") = sampens,
                            Rcpp::Named("meths") = meths);
}