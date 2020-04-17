#ifdef STAN_OPENCL
#include <stan/math/opencl/is_matrix_cl.hpp>
#include <gtest/gtest.h>
#include <type_traits>

TEST(type_trait, matrix_cl) {
  using stan::is_matrix_cl;
  using stan::math::matrix_cl;
  EXPECT_TRUE((is_matrix_cl<matrix_cl<double, -1, -1>>::value));
  EXPECT_FALSE((is_matrix_cl<double>::value));
}
#endif
