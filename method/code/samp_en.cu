#include <array>
#include <cstdint>
#include <cuda_runtime.h>
#include <iostream>
#include <vector>

#include "methData.h"
#include "samp_en.h"

__global__ void templateMatcher(int8_t *data, const unsigned int cumulativeSize,
                                const unsigned int m, unsigned int *d_prefixSum,
                                const unsigned int numBins, unsigned int *cm, unsigned int *cm_1) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i >= cumulativeSize - m || numBins == 0)
    return;
  unsigned int low  = 0;
  unsigned int high = numBins - 1;
  unsigned int mid;
  int          binIndex = -1;

  /// use binary search to find the bin to which the position belongs
  while (low <= high) {
    mid = low + (high - low) / 2;
    if (d_prefixSum[mid] == i) {
      binIndex = mid;
      break;
    } else if (d_prefixSum[mid] < i) {
      binIndex = mid;
      low      = mid + 1;
    } else {
      high = mid - 1;
    }
  }
  unsigned int binEnd  = (binIndex == numBins - 1) ? cumulativeSize : d_prefixSum[binIndex + 1];
  unsigned int binSize = binEnd - d_prefixSum[binIndex];
  unsigned int i_rel   = i - d_prefixSum[binIndex];

  if (i_rel >= binSize - m)
    return;

  bool         add = true;
  unsigned int idx = 0;
  for (unsigned int k = 0; k < m; k++) {
    if (data[i + k] == -1) {
      add = false;
      break;
    }
    idx += (data[i + k]) * (1 << k);
  }
  if (add && data[i + m] != -1) {
    atomicAdd(&cm[(1 << m) * binIndex + idx], 1);
    idx += (data[i + m]) * (1 << m);
    atomicAdd(&cm_1[(1 << (m + 1)) * binIndex + idx], 1);
  }
}

/**
 * Compute sample entropies at a bin level and aggregated at a file level for multiple files
 * simultaneously.
 *
 * @param fileMap FileMap object storing parsed methylation data of all files
 * @param m base length of templates
 * @return SampEns object storing per bin and file aggregated entropies
 */
SampEns sampEn(FileMap &fileMap, const unsigned int m, unsigned int n_streams) {
  SampEns sampens;
  Counts  counts;

  /**
   * initialise sampens with the default value of -1.0 which indicates that the sample entropy was
   * not computable
   * also initialise counts which keeps track of the total count of matching m-length and
   * (m+1)-length templates in every file
   */
  for (const auto &[file, fileMeths] : fileMap) {
    std::vector<std::vector<double>> x;
    for (const auto &chrBins : fileMeths)
      x.emplace_back(std::vector<double>(chrBins.meth.size(), -1.0));
    sampens[file] = FileSampEns{std::move(x), -1.0};
    counts[file]  = FileCounts{0, 0};
  }

  /// create CUDA streams
  cudaStream_t streams[n_streams];
  for (unsigned short i = 0; i < n_streams; ++i) {
    cudaStreamCreate(&streams[i]);
  }

  /**
   * array of pointers where each pointer is responsible for all methylation information for a
   * chromosome in a file
   */
  int8_t *d_flatBins[n_streams];
  /**
   * array of pointers where each pointer is responsible for tracking the number of m-length
   * templates of different types for a chromosome in a file
   */
  unsigned int *d_cm[n_streams];
  /**
   * array of pointers where each pointer is responsible for tracking the number of (m+1)-length
   * templates of different types for a chromosome in a file
   */
  unsigned int *d_cm_1[n_streams];
  /**
   * array of pointers where each pointer is responsible for tracking the starting positions of non
   * empty bins in d_flatBins
   */
  unsigned int *d_prefix_sum[n_streams];
  /**
   * array of pointers where each pointer is responsible for tracking the original indices of non
   * empty bins in a chromosome of a file
   */
  std::vector<unsigned int> goodBins[n_streams] = {};
  /**
   * array of pointers where each pointer is responsible for tracking the total number of positions
   * in a chromosome of a file
   */
  unsigned int cumulativeSize[n_streams] = {0};
  /**
   * array of pairs keeping track of the current filenames and chromosomes being worked on by the
   * different streams
   */
  std::pair<std::string, unsigned int> streamInfo[n_streams] = {{"", 0}};

  unsigned int   fileIndex = 0;
  unsigned short streamIdx = 0;

  for (const auto &[file, fileMeths] : fileMap) {
    for (unsigned int chrIndex = 0; chrIndex < fileMeths.size(); chrIndex++) {
      streamInfo[streamIdx] = {file, chrIndex};
      std::vector<unsigned int> prefix_sum;

      /// store the starting indices of good bins for a stream
      for (unsigned int i = 0; i < fileMeths[chrIndex].meth.size(); i++) {
        if (!fileMeths[chrIndex].meth[i].empty()) {
          goodBins[streamIdx].push_back(i);
          prefix_sum.push_back(cumulativeSize[streamIdx]);
          cumulativeSize[streamIdx] += fileMeths[chrIndex].meth[i].size();
        }
      }

      cudaMallocAsync(&d_flatBins[streamIdx], cumulativeSize[streamIdx] * sizeof(int8_t),
                      streams[streamIdx]);
      cudaMallocAsync(&d_prefix_sum[streamIdx], prefix_sum.size() * sizeof(unsigned int),
                      streams[streamIdx]);
      cudaMallocAsync(&d_cm[streamIdx],
                      goodBins[streamIdx].size() * (1 << m) * sizeof(unsigned int),
                      streams[streamIdx]);
      cudaMallocAsync(&d_cm_1[streamIdx],
                      goodBins[streamIdx].size() * (1 << (m + 1)) * sizeof(unsigned int),
                      streams[streamIdx]);

      /// store methylation information of all the bins in a flattened array
      for (unsigned int i = 0; i < goodBins[streamIdx].size(); i++) {
        unsigned int rowIndex = goodBins[streamIdx][i];
        cudaMemcpyAsync(d_flatBins[streamIdx] + prefix_sum[i],
                        fileMeths[chrIndex].meth[rowIndex].data(),
                        fileMeths[chrIndex].meth[rowIndex].size() * sizeof(int8_t),
                        cudaMemcpyHostToDevice, streams[streamIdx]);
      }

      cudaMemcpyAsync(d_prefix_sum[streamIdx], prefix_sum.data(),
                      prefix_sum.size() * sizeof(unsigned int), cudaMemcpyHostToDevice,
                      streams[streamIdx]);
      cudaMemsetAsync(d_cm[streamIdx], 0,
                      goodBins[streamIdx].size() * (1 << m) * sizeof(unsigned int),
                      streams[streamIdx]);
      cudaMemsetAsync(d_cm_1[streamIdx], 0,
                      goodBins[streamIdx].size() * (1 << (m + 1)) * sizeof(unsigned int),
                      streams[streamIdx]);

      unsigned int threadsPerBlock = 1024;
      unsigned int numBlocks = (cumulativeSize[streamIdx] + threadsPerBlock) / threadsPerBlock;

      templateMatcher<<<numBlocks, threadsPerBlock, 0, streams[streamIdx]>>>(
          d_flatBins[streamIdx], cumulativeSize[streamIdx], m, d_prefix_sum[streamIdx],
          goodBins[streamIdx].size(), d_cm[streamIdx], d_cm_1[streamIdx]);

      streamIdx++;

      if (streamIdx == n_streams ||
          (fileIndex == fileMap.size() - 1 && chrIndex == fileMeths.size() - 1)) {
        for (unsigned int j = 0; j < streamIdx; j++) {
          cudaStreamSynchronize(streams[j]);
          unsigned int *h_cm   = new unsigned int[goodBins[j].size() * (1 << m)]();
          unsigned int *h_cm_1 = new unsigned int[goodBins[j].size() * (1 << (m + 1))]();

          cudaMemcpy(h_cm, d_cm[j], goodBins[j].size() * (1 << m) * sizeof(unsigned int),
                     cudaMemcpyDeviceToHost);
          cudaMemcpy(h_cm_1, d_cm_1[j], goodBins[j].size() * (1 << (m + 1)) * sizeof(unsigned int),
                     cudaMemcpyDeviceToHost);

          // std::cout << "file: " << streamInfo[j].first << std::endl;
          // std::cout << "  chr: " << fileMap[streamInfo[j].first][streamInfo[j].second].chr <<
          // std::endl;

          for (unsigned int i = 0; i < goodBins[j].size(); i++) {
            // std::cout << "    Bin " << i << ":" << std::endl;
            // std::cout << "      cm:" << h_cm[i] << std::endl;
            // std::cout << "      cm_1:" << h_cm_1[i] << std::endl;
            unsigned long long cm   = 0;
            unsigned long long cm_1 = 0;

            /// compute the number of matching m-length and (m+1)-length templates
            for (unsigned int k = 0; k < (1 << (m + 1)); k++) {
              if (k < (1 << m) && h_cm[(1 << m) * i + k] > 1) {
                cm += ((unsigned long long)h_cm[(1 << m) * i + k] *
                       (unsigned long long)(h_cm[(1 << m) * i + k] - 1)) /
                      2;
              }
              if (h_cm_1[(1 << (m + 1)) * i + k] > 1) {
                cm_1 += ((unsigned long long)(h_cm_1[(1 << (m + 1)) * i + k]) *
                         (unsigned long long)(h_cm_1[(1 << (m + 1)) * i + k] - 1)) /
                        2;
              }
            }
            /**
             * compute sample entropy only when the number of matching m-length and (m+1)-length
             * templates are non-zero
             */
            if (cm != 0 && cm_1 != 0) {
              sampens[streamInfo[j].first].raw[streamInfo[j].second][goodBins[j][i]] =
                  log((double)cm / (double)cm_1);
              counts[streamInfo[j].first].cm += cm;
              counts[streamInfo[j].first].cm_1 += cm_1;
            }
          }
          // std::cout << std::endl
          //           << std::endl;

          cudaFreeAsync(d_flatBins[j], streams[j]);
          cudaFreeAsync(d_prefix_sum[j], streams[j]);
          cudaFreeAsync(d_cm[j], streams[j]);
          cudaFreeAsync(d_cm_1[j], streams[j]);
          goodBins[j].clear();
          cumulativeSize[j] = 0;

          delete[] h_cm;
          delete[] h_cm_1;
        }
        streamIdx = 0;
      }
    }
    fileIndex++;
  }
  // std::cout << std::endl;

  /// compute aggregated sample entropies at file level
  for (auto &[file, samp] : sampens) {
    if (counts[file].cm > 0 && counts[file].cm_1 > 0) {
      samp.agg = log((double)counts[file].cm / (double)counts[file].cm_1);
    }
  }

  return sampens;
}

// for (unsigned chrIndex = 0; chrIndex < meths.size(); chrIndex++)
// {
// std::vector<unsigned int> prefixSum;
// std::vector<unsigned int> goodBins;
// unsigned int *cm;
// unsigned int *cm_1;
// unsigned int *d_prefixSum;
// unsigned int cumulativeSize = 0;

// for (unsigned int i = 0; i < meths[chrIndex].meth.size(); i++)
// {
//   if (!meths[chrIndex].meth[i].empty())
//   {
//     goodBins.push_back(i);
//     prefixSum.push_back(cumulativeSize);
//     cumulativeSize += meths[chrIndex].meth[i].size();
//   }
// }

// cudaMalloc((void **)&d_flatBins, cumulativeSize * sizeof(char));

// for (unsigned int i = 0; i < goodBins.size(); i++)
// {
//   unsigned int rowIndex = goodBins[i];
//   cudaMemcpy(d_flatBins + prefixSum[i], meths[chrIndex].meth[rowIndex].data(),
//              meths[chrIndex].meth[rowIndex].size() * sizeof(char), cudaMemcpyHostToDevice);
// }

// cudaMalloc(&d_prefixSum, prefixSum.size() * sizeof(unsigned int));
// cudaMalloc((void **)&cm, goodBins.size() * sizeof(unsigned int));
// cudaMalloc((void **)&cm_1, goodBins.size() * sizeof(unsigned int));

// cudaMemset(cm, 0, goodBins.size() * sizeof(unsigned int));
// cudaMemset(cm_1, 0, goodBins.size() * sizeof(unsigned int));

// cudaMemcpy(d_prefixSum, prefixSum.data(), prefixSum.size() * sizeof(unsigned int),
// cudaMemcpyHostToDevice);

// dim3 threadsPerBlock(32, 32);
// dim3 numBlocks((cumulativeSize + threadsPerBlock.x) / threadsPerBlock.x, (cumulativeSize +
// threadsPerBlock.y) / threadsPerBlock.y);

// templateMatcher<<<numBlocks, threadsPerBlock>>>(d_flatBins, cumulativeSize, m, d_prefixSum,
// goodBins.size(), cm, cm_1);

// unsigned int *h_cm = new unsigned int[goodBins.size()];
// unsigned int *h_cm_1 = new unsigned int[goodBins.size()];

// cudaMemcpy(h_cm, cm, goodBins.size() * sizeof(unsigned int), cudaMemcpyDeviceToHost);
// cudaMemcpy(h_cm_1, cm_1, goodBins.size() * sizeof(unsigned int), cudaMemcpyDeviceToHost);

// cudaFree(d_flatBins);
// cudaFree(d_prefixSum);
// cudaFree(cm);
// cudaFree(cm_1);

// std::cout << "chr: " << meths[chrIndex].chr << std::endl;
// std::cout << "  cm: ";
// for (unsigned int i = 0; i < goodBins.size(); i++)
// {
//   std::cout << h_cm[i] << " ";
// }
// std::cout << std::endl;

// std::cout << "  cm_1: ";
// for (unsigned int i = 0; i < goodBins.size(); i++)
// {
//   std::cout << h_cm_1[i] << " ";
//   if (h_cm[i] != 0 && h_cm_1[i] != 0)
//     sampens[chrIndex][goodBins[i]] = log((double)h_cm[i] / (double)h_cm_1[i]);
// }
// std::cout << std::endl;
// }
// return sampens;

// DataRow *d_data;
// unsigned int *d_cm, *d_cm_1;
// unsigned int h_cm = 0, h_cm_1 = 0;

// cudaMalloc(&d_data, N * sizeof(DataRow));
// cudaMalloc(&d_cm, sizeof(unsigned int));
// cudaMalloc(&d_cm_1, sizeof(unsigned int));

// cudaMemcpy(d_data, data.data(), N * sizeof(DataRow), cudaMemcpyHostToDevice);
// cudaMemcpy(d_cm, &h_cm, sizeof(unsigned int), cudaMemcpyHostToDevice);
// cudaMemcpy(d_cm_1, &h_cm_1, sizeof(unsigned int), cudaMemcpyHostToDevice);

// dim3 threadsPerBlock(32, 32);
// dim3 numBlocks((N - m + threadsPerBlock.x) / threadsPerBlock.x, (N - m + threadsPerBlock.y) /
// threadsPerBlock.y);

// cudaEvent_t start, stop;
// cudaEventCreate(&start);
// cudaEventCreate(&stop);
// cudaEventRecord(start);
// templateMatcher<<<numBlocks, threadsPerBlock>>>(d_data, N, m, d_cm, d_cm_1);
// cudaEventRecord(stop);
// cudaEventSynchronize(stop);

// size_t free_mem, total_mem;
// cudaMemGetInfo(&free_mem, &total_mem);
// std::cout << "GPU Memory Used: " << (total_mem - free_mem) / (1024 * 1024) << "MB" <<
// std::endl;

// float milliseconds = 0;
// cudaEventElapsedTime(&milliseconds, start, stop);
// std::cout << "Sample Entropy computation time: " << milliseconds << " ms" << std::endl;

// cudaEventDestroy(start);
// cudaEventDestroy(stop);

// cudaMemcpy(&h_cm, d_cm, sizeof(unsigned int), cudaMemcpyDeviceToHost);
// cudaMemcpy(&h_cm_1, d_cm_1, sizeof(unsigned int), cudaMemcpyDeviceToHost);

// cudaFree(d_data);
// cudaFree(d_cm);
// cudaFree(d_cm_1);
// return log((double)h_cm / (double)h_cm_1);

// CPU only code below

// double sampEn(std::vector<DataRow> &data, const int m)
// {
//   unsigned int N = data.size();
//   unsigned int cm = 0, cm_1 = 0;
//   bool eq = true;
//   for (unsigned int i = 0; i < N - m; i++)
//   {
//     for (unsigned int j = i + 1; j < N - m; j++)
//     {
//       eq = true;
//       for (unsigned int k = 0; k < m; k++)
//       {
//         if (data[i + k].rate != data[j + k].rate)
//         {
//           eq = false;
//           break;
//         }
//       }
//       if (eq)
//         cm++;
//       if (eq && data[i + m].rate == data[j + m].rate)
//         cm_1++;
//     }
//   }
//   for (unsigned int i = 0; i < N - m; i++)
//   {
//     eq = true;
//     for (unsigned int k = 0; k < m; k++)
//     {
//       if (data[i + k].rate != data[N - m + k].rate)
//       {
//         eq = false;
//         break;
//       }
//     }
//     if (eq)
//       cm++;
//   }
//   return log((double)cm / (double)cm_1);
// }

// a transformation like
// unsigned int j = 1 + static_cast<unsigned int>((-1 + std::sqrt(8.0 * tid + 1)) / 2);
// unsigned int i = tid - (j * (j - 1)) / 2;
// could improve computational speed by considering the upper diagonal structure