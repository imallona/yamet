#include <iostream>
#include <vector>
#include <array>
#include <cuda_runtime.h>

#include "methData.h"
#include "samp_en.h"

__global__ void templateMatcher(char *data, const unsigned int cumulativeSize, const unsigned int m, unsigned int *d_prefixSum, const unsigned int numBins, unsigned int *cm, unsigned int *cm_1)
{
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int j = blockIdx.y * blockDim.y + threadIdx.y;

  if (j > cumulativeSize - m || i >= j || numBins == 0)
    return;
  unsigned int low = 0;
  unsigned int high = numBins - 1;
  unsigned int mid;
  int binIndex = -1;

  while (low <= high)
  {
    mid = low + (high - low) / 2;
    if (d_prefixSum[mid] == j)
    {
      binIndex = mid;
      break;
    }
    else if (d_prefixSum[mid] < j)
    {
      binIndex = mid;
      low = mid + 1;
    }
    else
    {
      high = mid - 1;
    }
  }
  unsigned int binEnd = (binIndex == numBins - 1) ? cumulativeSize : d_prefixSum[binIndex + 1];
  unsigned int binSize = binEnd - d_prefixSum[binIndex];
  unsigned int j_rel = j - d_prefixSum[binIndex];

  if (i < d_prefixSum[binIndex] || j_rel >= binSize - m)
    return;

  bool eq = true;
  for (unsigned int k = 0; k < m; k++)
  {
    if (data[i + k] == -1 || data[j + k] == -1 || data[i + k] != data[j + k])
    {
      eq = false;
      break;
    }
  }
  if (eq && data[i + m] != -1 && data[j + m] != -1)
  {
    atomicAdd(&cm[binIndex], 1);
    if (data[i + m] == data[j + m])
      atomicAdd(&cm_1[binIndex], 1);
  }
}

SampEns sampEn(FileMap &fileMap, const unsigned int m)
{
  SampEns sampens;
  Counts counts;
  for (const auto &[file, fileMeths] : fileMap)
  {
    std::vector<std::vector<double>> x;
    for (const auto &chrBins : fileMeths)
      x.emplace_back(std::vector<double>(chrBins.meth.size(), -1.0));
    sampens[file] = FileSampEns{std::move(x), -1.0};
    counts[file] = FileCounts{0, 0};
  }

  unsigned short num_streams = 50;

  // create CUDA streams
  cudaStream_t streams[num_streams];
  for (unsigned short i = 0; i < num_streams; ++i)
  {
    cudaStreamCreate(&streams[i]);
  }

  char *d_flatBins[num_streams];
  unsigned int *d_cm[num_streams];
  unsigned int *d_cm_1[num_streams];
  unsigned int *d_prefix_sum[num_streams];
  std::vector<unsigned int> goodBins[num_streams] = {};
  unsigned int cumulativeSize[num_streams] = {0};
  std::pair<std::string, unsigned int> streamInfo[num_streams] = {{"", 0}};

  unsigned int fileIndex = 0;
  unsigned short streamIdx = 0;

  for (const auto &[file, fileMeths] : fileMap)
  {
    for (unsigned int chrIndex = 0; chrIndex < fileMeths.size(); chrIndex++)
    {
      streamInfo[streamIdx] = {file, chrIndex};
      std::vector<unsigned int> prefix_sum;
      for (unsigned int i = 0; i < fileMeths[chrIndex].meth.size(); i++)
      {
        if (!fileMeths[chrIndex].meth[i].empty())
        {
          goodBins[streamIdx].push_back(i);
          prefix_sum.push_back(cumulativeSize[streamIdx]);
          cumulativeSize[streamIdx] += fileMeths[chrIndex].meth[i].size();
        }
      }

      cudaMallocAsync(&d_flatBins[streamIdx], cumulativeSize[streamIdx] * sizeof(char), streams[streamIdx]);
      cudaMallocAsync(&d_prefix_sum[streamIdx], prefix_sum.size() * sizeof(unsigned int), streams[streamIdx]);
      cudaMallocAsync(&d_cm[streamIdx], goodBins[streamIdx].size() * sizeof(unsigned int), streams[streamIdx]);
      cudaMallocAsync(&d_cm_1[streamIdx], goodBins[streamIdx].size() * sizeof(unsigned int), streams[streamIdx]);

      for (unsigned int i = 0; i < goodBins[streamIdx].size(); i++)
      {
        unsigned int rowIndex = goodBins[streamIdx][i];
        cudaMemcpyAsync(d_flatBins[streamIdx] + prefix_sum[i], fileMeths[chrIndex].meth[rowIndex].data(),
                        fileMeths[chrIndex].meth[rowIndex].size() * sizeof(char), cudaMemcpyHostToDevice, streams[streamIdx]);
      }

      cudaMemcpyAsync(d_prefix_sum[streamIdx], prefix_sum.data(), prefix_sum.size() * sizeof(unsigned int), cudaMemcpyHostToDevice, streams[streamIdx]);
      cudaMemsetAsync(d_cm[streamIdx], 0, goodBins[streamIdx].size() * sizeof(unsigned int), streams[streamIdx]);
      cudaMemsetAsync(d_cm_1[streamIdx], 0, goodBins[streamIdx].size() * sizeof(unsigned int), streams[streamIdx]);

      dim3 threadsPerBlock(32, 32);
      dim3 numBlocks((cumulativeSize[streamIdx] + threadsPerBlock.x) / threadsPerBlock.x, (cumulativeSize[streamIdx] + threadsPerBlock.y) / threadsPerBlock.y);

      templateMatcher<<<numBlocks, threadsPerBlock, 0, streams[streamIdx]>>>(d_flatBins[streamIdx], cumulativeSize[streamIdx], m, d_prefix_sum[streamIdx], goodBins[streamIdx].size(), d_cm[streamIdx], d_cm_1[streamIdx]);

      streamIdx++;

      if (streamIdx == num_streams || (fileIndex == fileMap.size() - 1 && chrIndex == fileMeths.size() - 1))
      {
        for (unsigned int j = 0; j < streamIdx; j++)
        {
          cudaStreamSynchronize(streams[j]);
          unsigned int *h_cm = new unsigned int[goodBins[j].size()]();
          unsigned int *h_cm_1 = new unsigned int[goodBins[j].size()]();

          cudaMemcpy(h_cm, d_cm[j], goodBins[j].size() * sizeof(unsigned int), cudaMemcpyDeviceToHost);
          cudaMemcpy(h_cm_1, d_cm_1[j], goodBins[j].size() * sizeof(unsigned int), cudaMemcpyDeviceToHost);

          // std::cout << "file: " << streamInfo[j].first << std::endl;
          // std::cout << "  chr: " << fileMap[streamInfo[j].first][streamInfo[j].second].chr << std::endl;

          for (unsigned int i = 0; i < goodBins[j].size(); i++)
          {
            // std::cout << "    Bin " << i << ":" << std::endl;
            // std::cout << "      cm:" << h_cm[i] << std::endl;
            // std::cout << "      cm_1:" << h_cm_1[i] << std::endl;
            if (h_cm[i] != 0 && h_cm_1[i] != 0)
            {
              sampens[streamInfo[j].first].raw[streamInfo[j].second][goodBins[j][i]] = log((double)h_cm[i] / (double)h_cm_1[i]);
              counts[streamInfo[j].first].cm += h_cm[i];
              counts[streamInfo[j].first].cm_1 += h_cm_1[i];
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

  for (auto &[file, samp] : sampens)
  {
    if (counts[file].cm > 0 && counts[file].cm_1 > 0)
    {
      samp.agg = log((double)counts[file].cm / (double)counts[file].cm_1);
    }
  }

  return sampens;

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

  // cudaMemcpy(d_prefixSum, prefixSum.data(), prefixSum.size() * sizeof(unsigned int), cudaMemcpyHostToDevice);

  // dim3 threadsPerBlock(32, 32);
  // dim3 numBlocks((cumulativeSize + threadsPerBlock.x) / threadsPerBlock.x, (cumulativeSize + threadsPerBlock.y) / threadsPerBlock.y);

  // templateMatcher<<<numBlocks, threadsPerBlock>>>(d_flatBins, cumulativeSize, m, d_prefixSum, goodBins.size(), cm, cm_1);

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
  // dim3 numBlocks((N - m + threadsPerBlock.x) / threadsPerBlock.x, (N - m + threadsPerBlock.y) / threadsPerBlock.y);

  // cudaEvent_t start, stop;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord(start);
  // templateMatcher<<<numBlocks, threadsPerBlock>>>(d_data, N, m, d_cm, d_cm_1);
  // cudaEventRecord(stop);
  // cudaEventSynchronize(stop);

  // size_t free_mem, total_mem;
  // cudaMemGetInfo(&free_mem, &total_mem);
  // std::cout << "GPU Memory Used: " << (total_mem - free_mem) / (1024 * 1024) << "MB" << std::endl;

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
}

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