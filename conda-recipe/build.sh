#!/bin/bash
set -euo pipefail

cd method

cmake -S. -Bbuild \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DVERSION="${PKG_VERSION}"

cmake --build build -j"${CPU_COUNT}" --config Release
cmake --install build
