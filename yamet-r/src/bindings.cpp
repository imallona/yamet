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
  Rcpp::NumericVector   cell_sampen(filenames.size());
  Rcpp::NumericVector   cell_meth(filenames.size());
  Rcpp::NumericVector   shannon(nrows);
  Rcpp::NumericVector   avg_meth(nrows);
  Rcpp::NumericMatrix   sampens(nrows, filenames.size());
  Rcpp::NumericMatrix   meths(nrows, filenames.size());

  unsigned int row = 0;
  for (unsigned int i = 0; i < intervals.size(); i++) {
    for (unsigned int j = 0; j < intervals[i].intervals.size(); j++) {
      unsigned int              total = 0;
      std::vector<unsigned int> shan((*fileMap.begin()).second.chrCounts[0].bins[0].cm.size(), 0);
      unsigned int              m_agg = 0, t_agg = 0;
      for (unsigned int k = 0; k < cpp_filenames.size(); k++) {
        sampens(row, k) = fileMap[cpp_filenames[k]].chrCounts[i].bins[j].sampen;
        if (sampens(row, k) == -1) {
          sampens(row, k) = NA_REAL;
        }
        meths(row, k) = fileMap[cpp_filenames[k]].chrCounts[i].bins[j].avg_meth;
        if (meths(row, k) == -1) {
          meths(row, k) = NA_REAL;
        }

        for (size_t l = 0; l < shan.size(); l++) {
          shan[l] += fileMap[cpp_filenames[k]].chrCounts[i].bins[j].cm[l];
          total += fileMap[cpp_filenames[k]].chrCounts[i].bins[j].cm[l];
        }
        m_agg += fileMap[cpp_filenames[k]].chrCounts[i].bins[j].m;
        t_agg += fileMap[cpp_filenames[k]].chrCounts[i].bins[j].t;
      }
      chr[row]   = intervals[i].chr;
      start[row] = intervals[i].intervals[j].start;
      end[row]   = intervals[i].intervals[j].end;

      double shannon_tmp = NA_REAL;
      if (total > 0) {
        shannon_tmp = 0;
        for (size_t l = 0; l < shan.size(); l++) {
          if (((double)shan[l]) / ((double)total) > 0) {
            shannon_tmp -=
                (((double)shan[l]) / ((double)total)) * log(((double)shan[l]) / ((double)total));
          }
        }
        shannon_tmp /= log(shan.size());
      }
      shannon[row] = shannon_tmp;

      double meth_tmp = NA_REAL;
      if (t_agg > 0) {
        meth_tmp = (((double)m_agg) / ((double)t_agg));
      }
      avg_meth[row] = meth_tmp;

      row++;
    }
  }

  for (unsigned int k = 0; k < cpp_filenames.size(); k++) {
    cell_sampen[k] = fileMap[cpp_filenames[k]].sampen;
    if (cell_sampen[k] == -1) {
      cell_sampen[k] = NA_REAL;
    }
    cell_meth[k] = fileMap[cpp_filenames[k]].avg_meth;
    if (cell_meth[k] == -1) {
      cell_sampen[k] = NA_REAL;
    }
  }

  return Rcpp::List::create(Rcpp::Named("chr") = chr, Rcpp::Named("start") = start,
                            Rcpp::Named("end") = end, Rcpp::Named("cell_sampen") = cell_sampen,
                            Rcpp::Named("cell_meth") = cell_meth, Rcpp::Named("shannon") = shannon,
                            Rcpp::Named("avg_meth") = avg_meth, Rcpp::Named("sampens") = sampens,
                            Rcpp::Named("meths") = meths);
}