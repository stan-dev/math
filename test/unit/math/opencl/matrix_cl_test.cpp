#ifdef STAN_OPENCL
#include <stan/math/mix/mat.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/fwd/mat/fun/typedefs.hpp>
#include <stan/math/mix/mat.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <vector>

TEST(MathMatrixCL, matrix_cl_double_creation) {
  stan::math::vector_d vec_1;
  stan::math::matrix_d mat_1;
  stan::math::matrix_d mat_2;

  vec_1.resize(3);
  mat_1.resize(2, 3);
  EXPECT_NO_THROW(stan::math::matrix_cl<double> A(1, 1));
  EXPECT_NO_THROW(stan::math::matrix_cl<double> vec_1cl(vec_1));
  EXPECT_NO_THROW(stan::math::matrix_cl<double> mat_1cl(mat_1));
  EXPECT_NO_THROW(stan::math::matrix_cl<double> mat_2cl(mat_2));
}


#endif
