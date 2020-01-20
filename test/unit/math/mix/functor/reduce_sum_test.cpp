#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

std::ostream* msgs = nullptr;

template <typename T>
T sum_(T arg) {
  return arg;
}

template <typename T, int RowType, int ColType>
auto sum_(const Eigen::Matrix<T, RowType, ColType>& arg) {
  return stan::math::sum(arg);
}

template <typename T>
auto sum_(const std::vector<T>& arg) {
  stan::scalar_type_t<T> sum = 0;
  for (size_t i = 0; i < arg.size(); ++i) {
    sum += sum_(arg[i]);
  }
  return sum;
}

struct sum_lpdf {
  template <typename T, typename... Args>
  inline auto operator()(std::size_t start, std::size_t end,
                         const std::vector<T>& sub_slice, std::ostream* msgs,
                         const Args&... args) const {
    using return_type = stan::return_type_t<T, Args...>;

    return stan::math::sum(sub_slice)
           + sub_slice.size()
                 * stan::math::sum(
                       std::vector<return_type>{return_type(sum_(args))...});
  }
};

TEST(MathMix_reduce_sum, double_slice) {
  auto f = [](const auto& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 0, msgs);
  };

  std::vector<double> data(5, 10.0);

  stan::test::expect_ad(f, data);
}

auto fi = [](const auto&... args) {
  return stan::math::reduce_sum<sum_lpdf>(std::vector<int>(2, 10.0), 0, msgs,
                                          args...);
};

auto fd = [](const auto& data, const auto&... args) {
  return stan::math::reduce_sum<sum_lpdf>(data, 0, msgs, args...);
};

TEST(MathMix_reduce_sum, int_arg) {
  std::vector<double> data(2, 1.0);
  int arg = 5.0;

  stan::test::expect_ad([&](const auto& data) { return fd(data, arg); }, data);
}

TEST(MathMix_reduce_sum, double_arg) {
  std::vector<double> data(2, 10.0);
  double arg = 5.0;

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_int_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<int> arg(2, 10);

  stan::test::expect_ad([&](const auto& data) { return fd(data, arg); }, data);
}

TEST(MathMix_reduce_sum, std_vector_double_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<double> arg(2, 10.0);

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(2);

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(2);

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, eigen_matrix_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_double_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<double>> arg(2, std::vector<double>(2, 10.0));

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(2));

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, eigen_three_args1) {
  Eigen::VectorXd arg1 = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd arg2 = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(fi, arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args2) {
  double arg1 = 1.0;
  std::vector<double> arg2(2, 1.0);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(fi, arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args3) {
  double arg1 = 1.0;
  std::vector<std::vector<double>> arg2(2, std::vector<double>(2, 1.0));
  std::vector<Eigen::MatrixXd> arg3(2, Eigen::MatrixXd::Ones(2, 2));

  stan::test::expect_ad(fi, arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_ints1) {
  Eigen::VectorXd arg1 = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd arg2 = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](const auto& arg1, const auto& arg2, const auto& arg3) {
        return fi(1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_ints2) {
  double arg1 = 1.0;
  std::vector<double> arg2(2, 1.0);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](const auto& arg1, const auto& arg2, const auto& arg3) {
        return fi(1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_ints3) {
  double arg1 = 1.0;
  std::vector<std::vector<double>> arg2(2, std::vector<double>(2, 1.0));
  std::vector<Eigen::MatrixXd> arg3(2, Eigen::MatrixXd::Ones(2, 2));

  stan::test::expect_ad(
      [&](const auto& arg1, const auto& arg2, const auto& arg3) {
        return fi(1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_doubles1) {
  Eigen::VectorXd arg1 = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd arg2 = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](const auto& arg1, const auto& arg2, const auto& arg3) {
        return fi(std::vector<double>{1.0, 2.0, 3.0}, arg1, 3.0, arg2,
                  std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_doubles2) {
  double arg1 = 1.0;
  std::vector<double> arg2(2, 1.0);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](const auto& arg1, const auto& arg2, const auto& arg3) {
        return fi(std::vector<double>{1.0, 2.0, 3.0}, arg1, 3.0, arg2,
                  std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_doubles3) {
  double arg1 = 1.0;
  std::vector<std::vector<double>> arg2(2, std::vector<double>(2, 1.0));
  std::vector<Eigen::MatrixXd> arg3(2, Eigen::MatrixXd::Ones(2, 2));

  stan::test::expect_ad(
      [&](const auto& arg1, const auto& arg2, const auto& arg3) {
        return fi(std::vector<double>{1.0, 2.0, 3.0}, arg1, 3.0, arg2,
                  std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
}
