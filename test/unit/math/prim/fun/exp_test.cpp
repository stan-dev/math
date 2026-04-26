

#include <stan/math.hpp>
#include <gtest/gtest.h>

#ifdef STAN_THREADS
TEST(MathFunctions, expInt) {
  using stan::math::exp;
  EXPECT_FLOAT_EQ(std::exp(3), exp(3));
  EXPECT_FLOAT_EQ(std::exp(3.1), exp(3.1));
  EXPECT_FLOAT_EQ(std::exp(3.0), exp(3.0));
}

TEST(MathFunctions, expVec) {
  using stan::math::exp;
  size_t N = 100;
  std::vector<double> vec(N);
  for (size_t i = 0; i < N; ++i) {
    vec[i] = i + 1;
  }
  std::vector<double> vec_test;
  EXPECT_NO_THROW(vec_test = stan::math::exp_test(vec));

  // std::vector<double> vec_test;
  // vec_test = stan::math::exp_test(vec);
  // for (size_t i = 0; i < N; ++i) {
  //   std::cout << vec_test[i] << "\n";
  //   std::cout << std::exp(i + 1) << "\n";
  //   //EXPECT_FLOAT_EQ(std::exp(i + 1), vec_test[i]);
  // }
}

#else
TEST(MathFunctions, expInt) {
  using stan::math::exp;
  EXPECT_FLOAT_EQ(std::exp(3), exp(3));
  EXPECT_FLOAT_EQ(std::exp(3.1), exp(3.1));
  EXPECT_FLOAT_EQ(std::exp(3.0), exp(3.0));
}
#endif
