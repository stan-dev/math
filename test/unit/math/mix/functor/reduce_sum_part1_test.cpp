#include <stan/math/prim/meta.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/functor/reduce_sum_util.hpp>

#include <limits>
#include <vector>

// Reduce sum tests are broken up into four files to avoid windows compiler
// error

TEST(MathMix_reduce_sum, grainsize_static) {
  using stan::math::test::get_new_msg;
  using stan::math::test::sum_lpdf;

  auto f1 = [](auto&& data) {
    return stan::math::reduce_sum_static<sum_lpdf>(data, 0, get_new_msg());
  };
  auto f2 = [](auto&& data) {
    return stan::math::reduce_sum_static<sum_lpdf>(data, -1, get_new_msg());
  };
  auto f3 = [](auto&& data) {
    return stan::math::reduce_sum_static<sum_lpdf>(data, 1, get_new_msg());
  };
  auto f4 = [](auto&& data) {
    return stan::math::reduce_sum_static<sum_lpdf>(data, 20, get_new_msg());
  };

  std::vector<double> data(5, 10.0);

  stan::test::expect_ad(f1, data);
  stan::test::expect_ad(f2, data);
  stan::test::expect_ad(f3, data);
  stan::test::expect_ad(f4, data);
}

TEST(MathMix_reduce_sum, grainsize) {
  using stan::math::test::get_new_msg;
  using stan::math::test::sum_lpdf;
  auto f1 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 0, get_new_msg());
  };
  auto f2 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, -1, get_new_msg());
  };
  auto f3 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, get_new_msg());
  };
  auto f4 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 20, get_new_msg());
  };

  std::vector<double> data(5, 10.0);

  stan::test::expect_ad(f1, data);
  stan::test::expect_ad(f2, data);
  stan::test::expect_ad(f3, data);
  stan::test::expect_ad(f4, data);
}

TEST(MathMix_reduce_sum, std_vector_zero_length) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<double> data(0);

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<double> data(5, 10.0);

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<std::vector<double>> data(3, std::vector<double>(2, 10.0));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_eigen_vector_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<Eigen::VectorXd> data(3, Eigen::VectorXd::Ones(2));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_eigen_row_vector_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<Eigen::RowVectorXd> data(3, Eigen::RowVectorXd::Ones(2));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_eigen_matrix_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<Eigen::MatrixXd> data(3, Eigen::MatrixXd::Ones(2, 4));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_std_vector_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<std::vector<std::vector<double>>> data(
      3, std::vector<std::vector<double>>(2, std::vector<double>(2, 10.0)));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_vector_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<std::vector<Eigen::VectorXd>> data(
      3, std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(2)));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_row_vector_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<std::vector<Eigen::RowVectorXd>> data(
      3, std::vector<Eigen::RowVectorXd>(2, Eigen::RowVectorXd::Ones(2)));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_matrix_double_slice) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;

  std::vector<std::vector<Eigen::MatrixXd>> data(
      3, std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Ones(2, 4)));

  stan::test::expect_ad(reduce_sum_static_sum_lpdf, data);
  stan::test::expect_ad(reduce_sum_sum_lpdf, data);
}

TEST(StanMath_reduce_sum_static, start_end_slice) {
  using stan::math::test::get_new_msg;
  using stan::math::test::start_end_lpdf;
  auto start_end_static = [](auto&& arg) {
    return stan::math::reduce_sum_static<start_end_lpdf>(arg, 1, get_new_msg(),
                                                         arg);
  };

  auto start_end = [](auto&& arg) {
    return stan::math::reduce_sum_static<start_end_lpdf>(arg, 1, get_new_msg(),
                                                         arg);
  };

  std::vector<double> data(5, 1.0);

  stan::test::expect_ad(start_end, data);
  stan::test::expect_ad(start_end_static, data);
}

TEST(MathMix_reduce_sum, int_arg) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;
  std::vector<double> data(2, 1.0);
  int arg = 5.0;

  stan::test::expect_ad(
      [&](auto&& data) { return reduce_sum_static_sum_lpdf(data, arg); }, data);
  stan::test::expect_ad(
      [&](auto&& data) { return reduce_sum_sum_lpdf(data, arg); }, data);
}

TEST(MathMix_reduce_sum, std_vector_int_arg) {
  using stan::math::test::reduce_sum_static_sum_lpdf;
  using stan::math::test::reduce_sum_sum_lpdf;
  std::vector<double> data(2, 10.0);
  std::vector<int> arg(2, 10);

  stan::test::expect_ad(
      [&](auto&& data) { return reduce_sum_static_sum_lpdf(data, arg); }, data);
  stan::test::expect_ad(
      [&](auto&& data) { return reduce_sum_sum_lpdf(data, arg); }, data);
}

TEST(MathMix_reduce_sum, double_arg) {
  stan::math::test::expect_ad_reduce_sum_lpdf(std::vector<double>(2, 10.0),
                                              5.0);
}

TEST(MathMix_reduce_sum, std_vector_double_arg) {
  stan::math::test::expect_ad_reduce_sum_lpdf(std::vector<double>(2, 10.0),
                                              std::vector<double>(2, 10.0));
}
