#ifndef STAN_MATH_OPENCL_OPENCL
#define STAN_MATH_OPENCL_OPENCL
#ifdef STAN_OPENCL

#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/add.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/copy_triangular.hpp>
#include <stan/math/opencl/cholesky_decompose.hpp>
#include <stan/math/opencl/diagonal_multiply.hpp>
#include <stan/math/opencl/identity.hpp>
#include <stan/math/opencl/tri_inverse.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/multiply_transpose.hpp>
#include <stan/math/opencl/sub_block.hpp>
#include <stan/math/opencl/subtract.hpp>
#include <stan/math/opencl/triangular_transpose.hpp>
#include <stan/math/opencl/transpose.hpp>
#include <stan/math/opencl/zeros.hpp>

#include <stan/math/opencl/err/check_diagonal_zeros.hpp>
#include <stan/math/opencl/err/check_matching_dims.hpp>
#include <stan/math/opencl/err/check_nan.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/opencl/err/check_square.hpp>
#include <stan/math/opencl/err/check_symmetric.hpp>
#endif
#endif
