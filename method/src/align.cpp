#include <atomic>
#include <cerrno>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#if defined(__linux__) || defined(__APPLE__)
#include <pthread.h>
// #include <sched.h>
#if defined(__APPLE__)
#include <mach/thread_act.h> // For pthread_mach_thread_np
// #include <mach/thread_policy.h> // For macOS-specific thread affinity API
#endif
#elif defined(_WIN32) // Windows
#include <windows.h>
#endif

#include <align.h>
#include <chrData.h>
#include <file_classes.h>
#include <window.h>

/**
 * Parse a tab separated file of all covered positions of a reference genome into a nested structure
 * from which the sample entropies can be calculated. The function is called in an asynchronous
 * manner for different files.
 *
 * @param filename[in] path to tab separated file to be parsed.
 * @param ref[in] Reference object with all positions of interest in the genome.
 * @param m[in] size of window.
 * @param chunk_size[in] size of file chunk.
 * @param fileMap[out] FileMap object of key value pairs where the parsed methylation data is stored
 * in a FileMeths object against the filename as key.
 */
void alignSingleWithRef(const std::string &filename, const Reference &ref, const unsigned int m,
                        const unsigned int chunk_size, FileMap &fileMap) {
  gzFile file = gzopen(filename.c_str(), "rb");
  if (!file) {
    throw std::system_error(errno, std::generic_category(), "Opening " + filename);
  }

  /**
   * buffer size of 64K by default - this is the total amount of information from a file stored at
   * any time as a chunk
   */
  std::vector<char> buffer(chunk_size);
  std::string       fullBuffer;
  // bool          headerSkipped = false;

  unsigned int chrIndex = 0, binIndex = 0, posIndex = 0;
  std::string  currentChr = "";
  bool         foundChr   = false;
  Window       window(m + 1);

  FileCounts fileCounts(ref, m);
  /**
   * for keeping track of the last covered position discovered, used so that we can register
   * breaks between covered positions, if any
   */
  int lastBinPos = -1;

  while (true) {
    int bytesRead = gzread(file, buffer.data(), chunk_size);

    if (bytesRead < 0) {
      gzclose(file);
      throw std::system_error(errno, std::generic_category(), filename);
    }

    fullBuffer.append(buffer.begin(), buffer.begin() + bytesRead);

    std::istringstream ss(fullBuffer);
    std::string        line;

    /// iterate over lines in a chunk
    while (std::getline(ss, line)) {
      /// save partial line to be processed later
      if (ss.peek() == EOF) {
        if (fullBuffer.back() == '\n') {
          fullBuffer.clear();
        } else {
          fullBuffer = line;
          break;
        }
      }
      // if (!headerSkipped) {
      //   headerSkipped = true;
      //   continue;
      // }
      std::istringstream lineStream(line);
      std::string        chr;
      unsigned int       pos;
      unsigned int       mVal, tVal, methValue;
      lineStream >> chr >> pos >> mVal >> tVal >> methValue;

      if (mVal > tVal) {
        throw std::system_error(
            EIO, std::generic_category(),
            "in line\n\n\t\033[33m" + line + "\033[0m\n\nof file " + filename +
                ". Total number of reads must exceed number of methylated reads.");
      }

      if (chr != currentChr) {
        binIndex = 0, posIndex = 0;
        window.clear();
        foundChr   = false;
        currentChr = chr;
      } else if (!foundChr) {
        continue;
      }

      if (!foundChr) {
        for (chrIndex = 0; chrIndex < ref.size(); chrIndex++) {
          if (ref[chrIndex].chr == chr) {
            foundChr = true;
            break;
          }
        }
        if (!foundChr) {
          continue;
        }
      }

      while (true) {
        bool increment = false, exit = false;
        if (ref[chrIndex].positions[binIndex].size() > 0) {
          if (pos == ref[chrIndex].positions[binIndex][posIndex]) {
            if (lastBinPos > -1 && posIndex - lastBinPos > 1) {
              window.append(-1).notify(fileCounts, chrIndex, binIndex);
            }
            window.append(methValue).notify(fileCounts, chrIndex, binIndex);
            fileCounts.addReads(mVal, tVal, chrIndex, binIndex);
            increment = true, exit = true;
            lastBinPos = posIndex;
          } else if (pos < ref[chrIndex].positions[binIndex][posIndex])
            exit = true;
          else {
            increment = true;
          }
        } else {
          if (binIndex == ref[chrIndex].positions.size() - 1) {
            exit = true;
          } else {
            binIndex += 1, posIndex = 0, lastBinPos = -1;
            window.clear();
          }
        }
        if (increment) {
          if (binIndex == ref[chrIndex].positions.size() - 1) {
            if (posIndex == ref[chrIndex].positions[binIndex].size() - 1) {
              lastBinPos = -1, exit = true;
            } else {
              posIndex++;
            }
          } else {
            if (posIndex == ref[chrIndex].positions[binIndex].size() - 1) {
              binIndex += 1, lastBinPos = -1, posIndex = 0;
              window.clear();
            } else {
              posIndex++;
            }
          }
        }
        if (exit) {
          break;
        }
      }
    }
    if (bytesRead < chunk_size) {
      break;
    }
  }
  gzclose(file);
  fileMap.addFile(filename, fileCounts.getContainer());
}

// Function to set thread affinity to a specific core
void set_thread_affinity(int core_id) {
#if defined(__linux__) // POSIX systems
  cpu_set_t cpuset;
  CPU_ZERO(&cpuset);
  CPU_SET(core_id, &cpuset);
  int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
  if (rc != 0) {
    std::cerr << "Error setting thread affinity: " << rc << "\n";
  }
#elif defined(__APPLE__)
  // macOS-specific implementation
  thread_affinity_policy policy = {core_id};
  thread_policy_set(pthread_mach_thread_np(pthread_self()), THREAD_AFFINITY_POLICY,
                    (thread_policy_t)&policy, 1);
#elif defined(_WIN32) // Windows
  DWORD_PTR mask = 1 << core_id;
  if (!SetThreadAffinityMask(GetCurrentThread(), mask)) {
    std::cerr << "Error setting thread affinity\n";
  }
#endif
}

/**
 * A wrapper function for alignSingleWithRef which allows multiple cell files to be parsed
 * simultaneously in a multi-threaded fashion
 *
 * @param filenames list of tab separated files to be parsed.
 * @param ref Reference object with all positions of interest in the genome.
 * @param m size of window.
 * @param chunk_size size of file chunk.
 * @return FileMap object of key value pairs where the parsed methylation data of every file is
 * stored in a FileMeths object against the filename as key
 */
FileMap alignWithRef(const std::vector<std::string> &filenames, const Reference &ref,
                     const unsigned int m, unsigned int n_cores, unsigned int n_threads_per_core,
                     const unsigned int chunk_size) {
  FileMap fileMap;
  fileMap.reserve(filenames.size());

  std::vector<std::thread> threads;
  std::atomic<bool>        done(false);
  std::atomic<bool>        exceptionOccurred(false);
  std::mutex               queueMutex;
  std::condition_variable  cv;
  std::exception_ptr       threadException;

  std::queue<std::string> taskQueue;
  for (const auto &filename : filenames) {
    taskQueue.push(filename);
  }

  for (unsigned i = 0; i < n_cores * n_threads_per_core; ++i) {
    threads.emplace_back([&, i] {
      try {
        // Set thread affinity
        set_thread_affinity(i / n_threads_per_core);

        // Process tasks
        while (true) {
          if (exceptionOccurred.load()) {
            return;
          }

          std::string filename;

          // Fetch task
          {
            std::unique_lock<std::mutex> lock(queueMutex);
            cv.wait(lock, [&]() { return done || !taskQueue.empty(); });

            if (taskQueue.empty()) {
              return; // Exit thread if no tasks remain
            }

            filename = taskQueue.front();
            taskQueue.pop();
          }
          // Process task
          alignSingleWithRef(filename, ref, m, chunk_size, fileMap);

          // Notify waiting threads
          cv.notify_all();
        }
      } catch (...) {
        exceptionOccurred.store(true);
        threadException = std::current_exception();
      }
    });
  }

  // Notify threads when all tasks are queued
  {
    std::unique_lock<std::mutex> lock(queueMutex);
    done = true;
  }
  cv.notify_all();

  // Join threads
  for (auto &t : threads) {
    if (t.joinable()) {
      t.join();
    }
  }

  if (exceptionOccurred.load() && threadException) {
    std::rethrow_exception(threadException);
  }

  return fileMap;
}