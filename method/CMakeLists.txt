# Minimum required CMake version
cmake_minimum_required(VERSION 3.15)

# Project name and version
project(yamet VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Boost Program Options and zlib
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(ZLIB REQUIRED)

# Add your source files (adjust the paths as needed)
file(GLOB SOURCES "code/*.cpp")
file(GLOB HEADERS "include/*.h")

# Define the executable target
add_executable(yamet ${SOURCES} ${HEADERS})

# Link Boost and zlib
target_link_libraries(yamet PUBLIC Boost::program_options ZLIB::ZLIB)

# Add include directories for headers
target_include_directories(yamet PUBLIC include)

# Install target (binary and headers)
install(TARGETS yamet DESTINATION bin)
install(DIRECTORY include/ DESTINATION include/yamet)
