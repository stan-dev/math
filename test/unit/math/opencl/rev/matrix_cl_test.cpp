#ifdef STAN_OPENCL
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/rev/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixRevCL, matrix_cl_var_creation) {
  using stan::math::matrix_cl;
  using stan::math::matrix_v;
  using stan::math::var;
  using stan::math::vector_v;
  vector_v vec_1;
  matrix_v mat_1;
  matrix_v mat_2;

  vec_1.resize(3);
  mat_1.resize(2, 3);
  mat_2.resize(3, 3);
  vec_1 << 1, 2, 3;
  mat_1 << 1, 2, 3, 1, 2, 3;
  mat_2 << 1, 2, 3, 1, 2, 3, 1, 2, 3;

  EXPECT_NO_THROW(matrix_cl<var> A(1, 1));
  EXPECT_NO_THROW(matrix_cl<var> d22(mat_1));
  EXPECT_NO_THROW(matrix_cl<var> d33(mat_2));
  EXPECT_NO_THROW(matrix_cl<var> d11(vec_1));
}

#endif
