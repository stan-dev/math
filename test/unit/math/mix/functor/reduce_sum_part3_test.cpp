#include <stan/math/prim/meta.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/functor/reduce_sum_util.hpp>

#include <limits>
#include <vector>

// Reduce sum tests are broken up into four files to avoid windows compiler
// error

TEST(MathMix_reduce_sum, eigen_three_args1) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  Eigen::VectorXd arg1 = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd arg2 = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(reduce_sum_static_int_sum_lpdf, arg1, arg2, arg3);
  stan::test::expect_ad(reduce_sum_int_sum_lpdf, arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args2) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  double arg1 = 1.0;
  std::vector<double> arg2(2, 1.0);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(reduce_sum_static_int_sum_lpdf, arg1, arg2, arg3);
  stan::test::expect_ad(reduce_sum_int_sum_lpdf, arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args3) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  double arg1 = 1.0;
  std::vector<std::vector<double>> arg2(2, std::vector<double>(2, 1.0));
  std::vector<Eigen::MatrixXd> arg3(2, Eigen::MatrixXd::Ones(2, 2));

  stan::test::expect_ad(reduce_sum_static_int_sum_lpdf, arg1, arg2, arg3);
  stan::test::expect_ad(reduce_sum_int_sum_lpdf, arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_ints1) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  Eigen::VectorXd arg1 = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd arg2 = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_static_int_sum_lpdf(
            1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_int_sum_lpdf(1, arg1, std::vector<int>{1, 2, 3}, arg2,
                                       3, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_ints2) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  double arg1 = 1.0;
  std::vector<double> arg2(2, 1.0);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_static_int_sum_lpdf(
            1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);
  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_int_sum_lpdf(1, arg1, std::vector<int>{1, 2, 3}, arg2,
                                       3, arg3);
      },
      arg1, arg2, arg3);
}
