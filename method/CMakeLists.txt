# Minimum required CMake version
cmake_minimum_required(VERSION 3.30)

# Project name and version
project(yamet VERSION 0.0.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Allow overriding the version from the command line
if(DEFINED VERSION)
  string(REGEX REPLACE "^v" "" VERSION "${VERSION}")
  set(PROJECT_VERSION ${VERSION})
  string(REGEX REPLACE "^([0-9]+).*" "\\1" PROJECT_VERSION_MAJOR ${VERSION})
  string(REGEX REPLACE "^[0-9]+\.([0-9]+).*" "\\1" PROJECT_VERSION_MINOR ${VERSION})
  string(REGEX REPLACE "^[0-9]+\.[0-9]+\.([0-9]+).*" "\\1" PROJECT_VERSION_PATCH ${VERSION})
endif()

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

# General package metadata
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Yet Another Methylation Entropy Tool")
set(CPACK_PACKAGE_VENDOR "Atreya Choudhury")
set(CPACK_PACKAGE_CONTACT "atreyachoudhury@hotmail.com")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/imallona/yamet")
set(CPACK_PACKAGE_LICENSE "GPL-3.0-only")

# Debian-specific metadata
set(CPACK_DEBIAN_PACKAGE_DEPENDS "zlib1g-dev (>= 1.2.0), libboost-program-options-dev (>= 1.70.0)")
set(CPACK_DEBIAN_PACKAGE_SECTION "utils")
set(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

# Optional: Specify the package output file name
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CMAKE_SYSTEM_NAME}-${CPACK_DEBIAN_ARCHITECTURE}-${VERSION}")

# Set version information
configure_file(
  ${CMAKE_SOURCE_DIR}/config/version.h.in
  ${CMAKE_BINARY_DIR}/config/version.h
)

# Add CPack configuration for packaging
include(CPack)