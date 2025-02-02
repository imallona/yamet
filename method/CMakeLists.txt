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

# Set up dependencies
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

# ============================
# Library Target Configuration
# ============================

# Library source file paths
file(GLOB LIB_SOURCES src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX "src/main.cpp|src/boost.cpp")

# Define the library target
add_library(${PROJECT_NAME}_lib STATIC ${LIB_SOURCES})

# Set library export details
set_target_properties(${PROJECT_NAME}_lib PROPERTIES
  OUTPUT_NAME ${PROJECT_NAME}
  EXPORT_NAME ${PROJECT_NAME}
)

# Link zlib and threads to library target
target_link_libraries(${PROJECT_NAME}_lib PRIVATE ZLIB::ZLIB Threads::Threads)

# Include directories for the library
target_include_directories(${PROJECT_NAME}_lib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/lib>
  $<INSTALL_INTERFACE:include>
)

# ===============================
# Executable Target Configuration
# ===============================

# Define the executable target
add_executable(${PROJECT_NAME} src/main.cpp src/boost.cpp)

# Link library and Boost program options to executable
target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_lib Boost::program_options)

# Link include directories for headers to executable target
target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include/cli>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/config>
)

# ============================
# Installation and Packaging
# ============================

# Include standard installation directories
include(GNUInstallDirs)

# Install library target
install(TARGETS ${PROJECT_NAME}_lib
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install executable target
install(TARGETS ${PROJECT_NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install public include folder for headers
install(DIRECTORY include/lib/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

# Export the library's targets
install(EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

include(CMakePackageConfigHelpers)

# Generate a config file for find_package
configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# Generate a config version file
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# Export targets file to build folder
export(EXPORT ${PROJECT_NAME}Targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
)

# Install the config and config version files
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}Config.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

# ===================
# CPack Configuration
# ===================

# General package metadata
set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Yet Another Methylation Entropy Tool")
set(CPACK_PACKAGE_VENDOR "Atreya Choudhury")
set(CPACK_PACKAGE_CONTACT "atreyachoudhury@hotmail.com")
set(CPACK_PACKAGE_HOMEPAGE_URL https://github.com/imallona/yamet)
set(CPACK_PACKAGE_LICENSE GPL-3.0-only)

# Debian-specific metadata
set(CPACK_DEBIAN_PACKAGE_DEPENDS "zlib1g-dev (>= 1.2.0)")
set(CPACK_DEBIAN_PACKAGE_SECTION "utils")
set(CPACK_DEBIAN_ARCHITECTURE ${CMAKE_SYSTEM_PROCESSOR})

# Package output file name
set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${CMAKE_SYSTEM_NAME}-${CPACK_DEBIAN_ARCHITECTURE}-${PROJECT_VERSION})

# Set version information
configure_file(
  ${CMAKE_SOURCE_DIR}/config/version.h.in
  ${CMAKE_BINARY_DIR}/config/version.h
)

# Add CPack configuration for packaging
include(CPack)