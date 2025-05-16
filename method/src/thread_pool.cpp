#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

using Task = std::function<void()>;

/**
 * This class manages a pool of worker threads to execute tasks concurrently. It provides a
 * mechanism to enqueue tasks and retrieve their results (or exceptions).
 */
class ThreadPool {
private:
  std::vector<std::thread>               workers;
  std::queue<std::packaged_task<void()>> taskQueue;
  std::mutex                             queueMutex;
  std::condition_variable                cv;
  std::atomic<bool>                      stop;
  std::atomic<bool>                      hasException;
  std::mutex                             exceptionMutex;
  std::exception_ptr                     storedException;

  /**
   * function executed by each worker thread.
   *
   * @param id The ID of the worker thread.
   */
  void worker_thread(int id) {
    while (!stop) {
      std::packaged_task<void()> task;
      {
        std::unique_lock<std::mutex> lock(queueMutex);
        cv.wait(lock, [this]() { return stop || !taskQueue.empty(); });
        if (stop && taskQueue.empty()) {
          return;
        }
        task = std::move(taskQueue.front());
        taskQueue.pop();
      }
      if (task.valid()) {
        try {
          task();
        } catch (...) {
          {
            std::lock_guard<std::mutex> lock(exceptionMutex);
            if (!storedException) {
              storedException = std::current_exception();
            }
            hasException.store(true);
          }
          cv.notify_all();
        }
      }
    }
  }

public:
  /**
   * Constructor for the ThreadPool.
   *
   * @param num_threads number of worker threads to create
   */
  explicit ThreadPool(size_t num_threads)
      : stop(false), hasException(false), storedException(nullptr) {
    for (size_t i = 0; i < num_threads; ++i) {
      workers.emplace_back(&ThreadPool::worker_thread, this, i);
    }
  }

  /**
   * Destructor for the ThreadPool.
   */
  ~ThreadPool() {
    stop = true;
    /// Signals all worker threads to stop and waits for them to finish.
    cv.notify_all();
    for (auto &worker : workers) {
      if (worker.joinable()) {
        worker.join();
      }
    }
  }

  /**
   * Enqueues a task for execution
   *
   * @param task function to be executed
   * @return std::future<void> that can be used to wait for the task
   * to complete and retrieve any potential exceptions
   */
  std::future<void> enqueue(Task task) {
    std::packaged_task<void()> packagedTask(std::move(task));
    std::future<void>          future = packagedTask.get_future();
    {
      std::unique_lock<std::mutex> lock(queueMutex);
      taskQueue.push(std::move(packagedTask));
    }
    cv.notify_one();
    return future;
  }

  /**
   * Rethrows the first exception caught by a worker thread, if any.
   *
   * @throws the stored exception if one exists
   */
  void rethrow_exception() {
    if (hasException.load() && storedException) {
      std::rethrow_exception(storedException);
    }
  }
};