#ifdef STAN_OPENCL
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/opencl.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixRevCL, matrix_cl_var_creation2) {
  using stan::math::matrix_cl;
  using stan::math::matrix_v;
  using stan::math::var;
  matrix_v mat_2;

  mat_2.resize(3, 3);
  mat_2 << 1, 2, 3, 1, 2, 3, 1, 2, 3;

  EXPECT_NO_THROW(matrix_cl<var> d33(mat_2));
}

#endif
