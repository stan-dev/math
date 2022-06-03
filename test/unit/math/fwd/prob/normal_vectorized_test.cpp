#include <stan/math/fwd.hpp>
#include <gtest/gtest.h>

TEST(ProbNormal, fvar_test_vlpdf) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<double>, -1, 1> Y = Eigen::Matrix<double, -1, 1>::Random(5);
  Eigen::Matrix<fvar<double>, -1, 1> Mu = Eigen::Matrix<double, -1, 1>::Random(5);
  Eigen::Matrix<fvar<double>, -1, 1> Sigma = stan::math::abs(Eigen::Matrix<double, -1, 1>::Random(5));
  Eigen::Matrix<fvar<double>, -1, 1> A = stan::math::normal_lpdf<false, stan::math::ProbReturnType::Vector>(Y, Mu, Sigma);

}
