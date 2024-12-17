# Minimum required CMake version
cmake_minimum_required(VERSION 3.15)

# Project name and version
project(yamet VERSION 0.1.0 LANGUAGES CXX)

# Check for the compiler and its version
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.0")
    message(FATAL_ERROR "GCC 8.0 or higher is required. Found: ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.0")
    message(FATAL_ERROR "Clang 8.0 or higher is required. Found: ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "8.0")
    message(FATAL_ERROR "AppleClang 8.0 or higher is required. Found: ${CMAKE_CXX_COMPILER_VERSION}")
  endif()
else()
  message(FATAL_ERROR "Unsupported compiler. This project requires GCC >= 8, Clang >= 8 or AppleClang >= 8")
endif()

# Set version information
configure_file(
  ${CMAKE_SOURCE_DIR}/config/version.h.in
  ${CMAKE_BINARY_DIR}/config/version.h
)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Boost Program Options and zlib
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(ZLIB REQUIRED)

# Add your source files (adjust the paths as needed)
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "include/*.h")

# Define the executable target
add_executable(yamet ${SOURCES} ${HEADERS})

# Link Boost and zlib
target_link_libraries(yamet PUBLIC Boost::program_options ZLIB::ZLIB)

# Add include directories for headers
target_include_directories(yamet PUBLIC include ${CMAKE_BINARY_DIR}/config)

# Install target (binary and headers)
install(TARGETS yamet DESTINATION bin)
install(DIRECTORY include/ DESTINATION include)