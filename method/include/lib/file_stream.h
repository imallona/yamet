#pragma once

#include <boost/iostreams/filtering_streambuf.hpp>

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
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

  CompressionType detectCompression(const std::string &filename);

public:
  explicit FileStream(const std::string &filename, const unsigned int chunk_size);
  ~FileStream();
  bool            good() const;
  void            close();
  std::streamsize read(char *buffer, std::streamsize size);
  bool            getline(std::string &line);
};