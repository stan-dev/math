#ifdef STAN_OPENCL
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingOpenCL, checkThrows) {
  const char* function = "test_func";
  const char* msg = "test";
  EXPECT_THROW(stan::math::throw_openCL(function, msg), std::domain_error);
}
#else
#include <gtest/gtest.h>
TEST(ErrorHandlingOpenCL, checkThrowsDummy) { EXPECT_NO_THROW(); }
#endif
