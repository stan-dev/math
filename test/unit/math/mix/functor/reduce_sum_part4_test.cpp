#include <stan/math/prim/meta.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/functor/reduce_sum_util.hpp>

#include <limits>
#include <vector>

// Reduce sum tests are broken up into four files to avoid windows compiler
// error

TEST(MathMix_reduce_sum, eigen_three_args_with_ints3) {
  using stan::math::test::reduce_sum_int_sum_lpdf;
  using stan::math::test::reduce_sum_static_int_sum_lpdf;
  double arg1 = 1.0;
  std::vector<std::vector<double>> arg2(2, std::vector<double>(2, 1.0));
  std::vector<Eigen::MatrixXd> arg3(2, Eigen::MatrixXd::Ones(2, 2));

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

TEST(MathMix_reduce_sum, eigen_three_args_with_doubles1) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;
  Eigen::VectorXd arg1 = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd arg2 = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_static_sum_lpdf(
            std::vector<double>{1.0, 2.0, 3.0}, arg1, 3.0, arg2,
            std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_sum_lpdf(std::vector<double>{1.0, 2.0, 3.0}, arg1,
                                   3.0, arg2,
                                   std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_doubles2) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;
  double arg1 = 1.0;
  std::vector<double> arg2(2, 1.0);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_static_sum_lpdf(
            std::vector<double>{1.0, 2.0, 3.0}, arg1, 3.0, arg2,
            std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_sum_lpdf(std::vector<double>{1.0, 2.0, 3.0}, arg1,
                                   3.0, arg2,
                                   std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_doubles3) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;
  double arg1 = 1.0;
  std::vector<std::vector<double>> arg2(2, std::vector<double>(2, 1.0));
  std::vector<Eigen::MatrixXd> arg3(2, Eigen::MatrixXd::Ones(2, 2));

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_static_sum_lpdf(
            std::vector<double>{1.0, 2.0, 3.0}, arg1, 3.0, arg2,
            std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return reduce_sum_sum_lpdf(std::vector<double>{1.0, 2.0, 3.0}, arg1,
                                   3.0, arg2,
                                   std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
}

#ifdef STAN_THREADS
TEST(MathMix_reduce_sum, static_check) {
  using stan::math::test::get_new_msg;
  using stan::math::test::static_check_lpdf;

  for (auto size : {1, 3, 6, 11}) {
    std::vector<int> data(size, 10);
    std::vector<double> arg(size, 10);

    auto fi1 = [&](auto&&... args) {
      return stan::math::reduce_sum_static<static_check_lpdf<1>>(
          data, 1, get_new_msg(), args...);
    };

    auto fi2 = [&](auto&&... args) {
      return stan::math::reduce_sum_static<static_check_lpdf<2>>(
          data, 2, get_new_msg(), args...);
    };

    auto fi3 = [&](auto&&... args) {
      return stan::math::reduce_sum_static<static_check_lpdf<3>>(
          data, 3, get_new_msg(), args...);
    };

    stan::test::expect_ad(fi1, arg);
    stan::test::expect_ad(fi2, arg);
    stan::test::expect_ad(fi3, arg);
  }
}
#endif
