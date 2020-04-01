#include <stan/math/prim/meta.hpp>
#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/reduce_sum_util.hpp>

#include <limits>
#include <vector>

TEST(MathMix_reduce_sum, grainsize) {
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

TEST(MathMix_reduce_sum, eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::VectorXd arg = Eigen::VectorXd::Ones(2);
  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::RowVectorXd arg = Eigen::RowVectorXd::Ones(2);

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, eigen_matrix_arg) {
  std::vector<double> data(2, 10.0);
  Eigen::MatrixXd arg = Eigen::MatrixXd::Ones(2, 2);

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_double_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<double>> arg(2, std::vector<double>(2, 10.0));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::VectorXd> arg(2, Eigen::VectorXd::Ones(2));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::RowVectorXd> arg(2, Eigen::RowVectorXd::Ones(2));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_eigen_matrix_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<Eigen::MatrixXd> arg(2, Eigen::MatrixXd::Ones(2, 2));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_std_vector_double_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<std::vector<double>>> arg(
      2, std::vector<std::vector<double>>(2, std::vector<double>(2, 10.0)));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<Eigen::VectorXd>> arg(
      2, std::vector<Eigen::VectorXd>(2, Eigen::VectorXd::Ones(2)));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_row_vector_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<Eigen::RowVectorXd>> arg(
      2, std::vector<Eigen::RowVectorXd>(2, Eigen::RowVectorXd::Ones(2)));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

TEST(MathMix_reduce_sum, std_vector_std_vector_eigen_matrix_arg) {
  std::vector<double> data(2, 10.0);
  std::vector<std::vector<Eigen::MatrixXd>> arg(
      2, std::vector<Eigen::MatrixXd>(2, Eigen::MatrixXd::Ones(2, 2)));

  stan::math::test::expect_ad_reduce_sum_lpdf(data, arg);
}

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

TEST(MathMix_reduce_sum, static_check) {
  stan::math::init_threadpool_tbb();
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
