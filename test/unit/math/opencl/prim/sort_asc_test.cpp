#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>


TEST(OpenCLPrim, sort_asc_simple) {
  Eigen::VectorXd a(9);
  a << 1, 2, 7, 4, 5, 6, 8, 3, 9;

  stan::math::matrix_cl<double> a_cl(a);

  Eigen::VectorXd correct = stan::math::sort_asc(a);
  Eigen::VectorXd res = stan::math::from_matrix_cl(stan::math::sort_asc(a_cl));

//  std::cout << res << std::endl << std::endl;
//  std::cout << correct << std::endl << std::endl;

  EXPECT_MATRIX_EQ(res, correct);
}

TEST(OpenCLPrim, sort_asc_large) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(37985);

  stan::math::matrix_cl<double> a_cl(a);

  Eigen::VectorXd correct = stan::math::sort_asc(a);
  Eigen::VectorXd res = stan::math::from_matrix_cl(stan::math::sort_asc(a_cl));

  EXPECT_MATRIX_EQ(res, correct);
}

#endif
