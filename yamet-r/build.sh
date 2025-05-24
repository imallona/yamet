#!/bin/bash

# Get the absolute path to the directory where the script is located
SCRIPT_DIR="../method"

# Set the install prefix (you can modify this to another path if needed)
INSTALL_PREFIX="./yamet"

# Set the build directory (relative to the current working directory)
BUILD_DIR="build"

# Create the build directory if it doesn't exist
if [ ! -d "$BUILD_DIR" ]; then
  mkdir "$BUILD_DIR"
fi

# Run CMake from the parent directory, pointing to the build directory
cmake -S"$SCRIPT_DIR" -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
  -DCMAKE_BUILD_TYPE=Release \
  -B"$BUILD_DIR"

# Build the project
cmake --build build --config Release

# Install the project to the specified location
cmake --install build

mkdir -p inst/include/yamet
mkdir -p inst/libs
cp -r $INSTALL_PREFIX/include/yamet/* inst/include/yamet
cp -r $INSTALL_PREFIX/lib/*yamet* inst/libs

Rscript -e 'Rcpp::compileAttributes(".")'
R CMD build --resave-data --compact-vignettes="both" .
