#!/bin/bash

YAMET_DIR="../method"

# Absolute path to the directory where the script is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Set the install prefix (you can modify this to another path if needed)
INSTALL_PREFIX="$SCRIPT_DIR/yamet"

# Set the build directory (relative to the current working directory)
BUILD_DIR="$SCRIPT_DIR/build"

# Create the build directory if it doesn't exist
if [ ! -d "$BUILD_DIR" ]; then
  mkdir "$BUILD_DIR"
fi

# Run CMake from the parent directory, pointing to the build directory
cmake -S"$YAMET_DIR" -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
  -DCMAKE_BUILD_TYPE=Release \
  -B"$BUILD_DIR"

# Build the project
cmake --build $BUILD_DIR -j 2 --config Release

# Install the project to the specified location
cmake --install $BUILD_DIR

mkdir -p $SCRIPT_DIR/inst/include/yamet
mkdir -p $SCRIPT_DIR/inst/libs
cp -r $INSTALL_PREFIX/include/yamet/* $SCRIPT_DIR/inst/include/yamet
cp -r $INSTALL_PREFIX/lib/*yamet* $SCRIPT_DIR/inst/libs

Rscript -e "Rcpp::compileAttributes('$SCRIPT_DIR')"
R CMD build --resave-data --compact-vignettes="both" $SCRIPT_DIR
