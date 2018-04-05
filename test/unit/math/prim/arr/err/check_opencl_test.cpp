#ifdef STAN_OPENCL
#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>

TEST(ErrorHandlingOpenCL, checkThrows) {
  const char* function = "test_func";
  cl::Error e(-5, "CL_OUT_OF_RESOURCES");
  EXPECT_THROW(stan::math::check_opencl_error(function, e), std::system_error);
}
#else
#include <gtest/gtest.h>
TEST(ErrorHandlingOpenCL, checkThrowsDummy) { EXPECT_NO_THROW(); }
#endif
