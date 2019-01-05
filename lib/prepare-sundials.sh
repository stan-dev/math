#!/bin/bash

SUNDIALS_VERSION="4.0.1"

SUNDIALS_TAR=sundials-${SUNDIALS_VERSION}.tar.gz

tar xvzf $SUNDIALS_TAR
mv sundials-${SUNDIALS_VERSION} sundials_${SUNDIALS_VERSION}

cd sundials_${SUNDIALS_VERSION}

rm -rf INSTALL_GUIDE.pdf config/ doc/ examples/ test/ */cvode */ida */arkode */kinsol
find . -name CMakeLists.txt -exec rm {} \;

find src -name "*.c" -type f -exec sed -E -i _orig  's#[^sf]printf\(#STAN_SUNDIALS_PRINTF(#'g {} \;
find src -name "*.c" -type f -exec sed -E -i _orig  's#fprintf\(#STAN_SUNDIALS_FPRINTF(#'g {} \;

find src -name "*.c_orig" -exec rm {} \;

echo "Please setup sundials_${SUNDIALS_VERSION}/include/sundials/sundials_config.h"
echo "Plese also ensure that Stan overrides are in place: sundials_${SUNDIALS_VERSION}/include/stan_sundials_printf_override.hpp"
