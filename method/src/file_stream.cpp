#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <cerrno>
#include <fstream>
#include <ios>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace io = boost::iostreams;

/**
 * Provides a stream interface for reading compressed files with automatic decompression.
 *
 * This class supports both Gzip (.gz) and Zstandard (.zst, .zstd) compressed files, as well as
 * uncompressed files. The compression type is automatically detected either by file extension or by
 * examining the file's magic numbers.
 */
class FileStream {
private:
  std::unique_ptr<std::ifstream>            file_;
  std::unique_ptr<io::filtering_istreambuf> filtering_buf_;
  std::unique_ptr<std::istream>             stream_;

  std::vector<char>  buffer_;
  std::string        fullLineBuffer_;
  std::istringstream currentLineStream_;
  unsigned int       chunkSize_;
  bool               eofReached_;

  enum class CompressionType { NONE, GZIP, ZSTD };

  /**
   * Detects the compression type of a file.
   * First checks the file extension, then falls back to examining the file's magic numbers if the
   * extension is not recognized.
   *
   * @param filename The name of the file to examine
   * @return The detected compression type
   */
  CompressionType detectCompression(const std::string &filename) {
    if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz") {
      return CompressionType::GZIP;
    }
    if (filename.size() >= 4 && filename.substr(filename.size() - 4) == ".zst") {
      return CompressionType::ZSTD;
    }
    if (filename.size() >= 5 && filename.substr(filename.size() - 5) == ".zstd") {
      return CompressionType::ZSTD;
    }

    std::ifstream test_file(filename, std::ios::binary);
    if (!test_file) {
      return CompressionType::NONE;
    }

    unsigned char header[4];
    test_file.read(reinterpret_cast<char *>(header), 4);
    test_file.close();

    /// gzip magic bytes: 1f 8b
    if (header[0] == 0x1f && header[1] == 0x8b) {
      return CompressionType::GZIP;
    }

    /// zstd magic bytes: 28 b5 2f fd
    if (header[0] == 0x28 && header[1] == 0xb5 && header[2] == 0x2f && header[3] == 0xfd) {
      return CompressionType::ZSTD;
    }

    return CompressionType::NONE;
  }

public:
  /**
   * Constructs a FileStream and opens the specified file
   * Automatically detects and configures the appropriate decompression strategy based on the file's
   * CompressionType obtained from `detectCompression`.
   *
   * @param filename The name of the file to open
   * @throw std::system_error if the file cannot be opened
   * @see detectCompression() For details on compression detection logic
   */
  explicit FileStream(const std::string &filename, const unsigned int chunk_size)
      : file_(std::make_unique<std::ifstream>(filename, std::ios::binary)), chunkSize_(chunk_size),
        buffer_(chunk_size), eofReached_(false) {
    if (!*file_) {
      throw std::system_error(errno, std::generic_category(), "Opening file " + filename);
    }

    filtering_buf_ = std::make_unique<io::filtering_istreambuf>();

    CompressionType compression = detectCompression(filename);

    switch (compression) {
    case CompressionType::GZIP:
      filtering_buf_->push(io::gzip_decompressor());
      break;
    case CompressionType::ZSTD:
      filtering_buf_->push(io::zstd_decompressor());
      break;
    case CompressionType::NONE:
      break;
    }

    filtering_buf_->push(*file_);
    stream_ = std::make_unique<std::istream>(filtering_buf_.get());
  }

  /**
   * Automatically closes the file if it's still open
   */
  ~FileStream() {
    close();
  }

  /**
   * Checks if the stream is in a good state
   * @return true if the stream is good, false otherwise
   */
  bool good() const {
    return stream_ && stream_->good();
  }

  /**
   * Closes the file and releases all resources
   */
  void close() {
    if (stream_) {
      stream_.reset();
    }
    if (filtering_buf_) {
      filtering_buf_.reset();
    }
    if (file_ && file_->is_open()) {
      file_->close();
    }
  }

  /**
   * Reads data from the decompressed stream
   *
   * @param buffer Pointer to the buffer where data will be stored
   * @param size Maximum number of bytes to read
   * @return The number of bytes actually read
   */
  std::streamsize read(char *buffer, std::streamsize size) {
    stream_->read(buffer, size);
    return stream_->gcount();
  }

  /**
   * Read the next line from the compressed file
   * @param line Reference to string where the line will be stored
   * @return true if a line was successfully read, false if EOF reached
   */
  bool getline(std::string &line) {
    line.clear();

    while (true) {
      if (std::getline(currentLineStream_, line)) {
        if (currentLineStream_.peek() == EOF) {
          currentLineStream_.clear();
          if (!fullLineBuffer_.empty() && fullLineBuffer_.back() == '\n') {
            fullLineBuffer_.clear();
            return true;
          } else {
            fullLineBuffer_ = line;
          }
        } else {
          return true;
        }
      } else {
        if (eofReached_) {
          if (!fullLineBuffer_.empty()) {
            line = fullLineBuffer_;
            fullLineBuffer_.clear();
            currentLineStream_.clear();
            return true;
          }
          return false;
        }

        std::streamsize bytesRead = read(buffer_.data(), chunkSize_);

        if (bytesRead < 0) {
          throw std::system_error(errno, std::generic_category(), "Error reading compressed file");
        }

        if (bytesRead == 0) {
          eofReached_ = true;
          continue;
        }

        fullLineBuffer_.append(buffer_.begin(), buffer_.begin() + bytesRead);

        currentLineStream_ = std::istringstream(fullLineBuffer_);

        if (bytesRead < chunkSize_) {
          eofReached_ = true;
        }
      }
    }
  }
};