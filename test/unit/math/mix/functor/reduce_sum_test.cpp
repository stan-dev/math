#include <stan/math/prim/meta.hpp>
#include <test/unit/math/test_ad.hpp>

#include <limits>
#include <vector>

std::ostream* msgs = nullptr;

template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
T sum_(T arg) {
  return arg;
}

template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
auto sum_(EigMat&& arg) {
  return stan::math::sum(arg);
}

template <typename Vec, stan::require_std_vector_t<Vec>* = nullptr>
auto sum_(Vec&& arg) {
  stan::scalar_type_t<Vec> sum = 0;
  for (size_t i = 0; i < arg.size(); ++i) {
    sum += sum_(arg[i]);
  }
  return sum;
}

struct sum_lpdf {
  template <typename T, typename... Args>
  inline auto operator()(std::size_t start, std::size_t end, T&& sub_slice,
                         std::ostream* msgs, Args&&... args) const {
    using return_type = stan::return_type_t<T, Args...>;

    return sum_(sub_slice)
           + sub_slice.size()
                 * stan::math::sum(std::vector<return_type>{
                       return_type(sum_(std::forward<Args>(args)))...});
  }
};

TEST(MathMix_reduce_sum, grainsize) {
  auto f1 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 0, msgs);
  };
  auto f2 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, -1, msgs);
  };
  auto f3 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };
  auto f4 = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 20, msgs);
  };

  std::vector<double> data(5, 10.0);

  stan::test::expect_ad(f1, data);
  stan::test::expect_ad(f2, data);
  stan::test::expect_ad(f3, data);
  stan::test::expect_ad(f4, data);
}

TEST(MathMix_reduce_sum, std_vector_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<double> data(5, 10.0);

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<std::vector<double>> data(3, std::vector<double>(2, 10.0));

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_eigen_vector_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<Eigen::VectorXd> data(3, Eigen::VectorXd::Ones(2));

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_eigen_row_vector_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<Eigen::RowVectorXd> data(3, Eigen::RowVectorXd::Ones(2));

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_eigen_matrix_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<Eigen::MatrixXd> data(3, Eigen::MatrixXd::Ones(2, 4));

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_std_vector_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<std::vector<std::vector<double>>> data(3, std::vector<std::vector<double>>(2, std::vector<double>(2, 10.0)));

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_vector_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<std::vector<Eigen::VectorXd>> data(3, std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(2)));

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_row_vector_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<std::vector<Eigen::RowVectorXd>> data(3, std::vector<Eigen::RowVectorXd>(2, Eigen::RowVectorXd::Ones(2)));

  stan::test::expect_ad(f, data);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_matrix_double_slice) {
  auto f = [](auto&& data) {
    return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs);
  };

  std::vector<std::vector<Eigen::MatrixXd>> data(3, std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Ones(2, 4)));

  stan::test::expect_ad(f, data);
}

struct start_end_lpdf {
  template <typename T1, typename T2>
  inline auto operator()(std::size_t start, std::size_t end, T1&&,
                         std::ostream* msgs, T2&& data) const {
    stan::return_type_t<T1, T2> sum = 0;
    EXPECT_GE(start, 0);
    EXPECT_LE(end, data.size() - 1);
    for (size_t i = start; i <= end; i++) {
      sum += data[i];
    }
    return sum;
  }
};

TEST(StanMath_reduce_sum, start_end_slice) {
  auto start_end = [](auto&& arg) {
    return stan::math::reduce_sum<start_end_lpdf>(arg, 1, msgs, arg);
  };

  std::vector<double> data(5, 1.0);

  stan::test::expect_ad(start_end, data);
}

auto fi = [](auto&&... args) {
  return stan::math::reduce_sum<sum_lpdf>(std::vector<int>(2, 10.0), 1, msgs,
                                          args...);
};

auto fd = [](auto&& data, auto&&... args) {
  return stan::math::reduce_sum<sum_lpdf>(data, 1, msgs, args...);
};

TEST(MathMix_reduce_sum, int_arg) {
  std::vector<double> data(2, 1.0);
  int arg = 5.0;

  stan::test::expect_ad([&](auto&& data) { return fd(data, arg); }, data);
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

  stan::test::expect_ad([&](auto&& data) { return fd(data, arg); }, data);
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

TEST(MathMix_reduce_sum, std_vector_eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::RowVectorXd> arg(2, Eigen::RowVectorXd::Ones(2));

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_eigen_matrix_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::MatrixXd> arg(2, Eigen::MatrixXd::Ones(2, 2));

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_std_vector_double_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<std::vector<double>>> arg(2, std::vector<std::vector<double>>(2, std::vector<double>(2, 10.0)));

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<Eigen::VectorXd>> arg(2, std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(2)));

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<Eigen::RowVectorXd>> arg(2, std::vector<Eigen::RowVectorXd>(2, Eigen::RowVectorXd::Ones(2)));

  stan::test::expect_ad(fi, arg);
  stan::test::expect_ad(fd, data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_matrix_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<Eigen::MatrixXd>> arg(2, std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Ones(2, 2)));

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
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return fi(1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_ints2) {
  double arg1 = 1.0;
  std::vector<double> arg2(2, 1.0);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return fi(1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_ints3) {
  double arg1 = 1.0;
  std::vector<std::vector<double>> arg2(2, std::vector<double>(2, 1.0));
  std::vector<Eigen::MatrixXd> arg3(2, Eigen::MatrixXd::Ones(2, 2));

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return fi(1, arg1, std::vector<int>{1, 2, 3}, arg2, 3, arg3);
      },
      arg1, arg2, arg3);
}

TEST(MathMix_reduce_sum, eigen_three_args_with_doubles1) {
  Eigen::VectorXd arg1 = Eigen::VectorXd::Ones(2);
  Eigen::RowVectorXd arg2 = Eigen::RowVectorXd::Ones(2);
  Eigen::MatrixXd arg3 = Eigen::MatrixXd::Ones(2, 2);

  stan::test::expect_ad(
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
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
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
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
      [&](auto&& arg1, auto&& arg2, auto&& arg3) {
        return fi(std::vector<double>{1.0, 2.0, 3.0}, arg1, 3.0, arg2,
                  std::vector<double>{1.0, 2.0, 3.0}, arg3);
      },
      arg1, arg2, arg3);
}
