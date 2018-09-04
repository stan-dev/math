#include <test/unit/math/rev/mat/util/test_autodiff.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <vector>

TEST(AgradRevMatrix, test_autodiff_scalar) {
  auto func = [](auto a, auto b) { return a + b; };

  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, 1, 1.0));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, 1.0, 1.0));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, 1.0, 1));
}

TEST(AgradRevMatrix, test_autodiff_std_vector) {
  auto func = [](auto a, auto b) { return stan::math::append_array(a, b); };

  std::vector<int> v1 = {{1, 2, 3}};
  std::vector<double> v2 = {{1.0, 2.0}};

  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v2));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v2));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v1));
}

TEST(AgradRevMatrix, test_autodiff_Eigen_VectorXd) {
  auto func = [](auto a, auto b) {
    stan::math::check_size_match("test_autodiff_Eigen_VectorXd", "a", a.size(),
                                 "b", b.size());
    return (a + b).eval();
  };

  Eigen::VectorXd v1(3);
  Eigen::VectorXd v2(2);

  v1 << 1.0, 2.0, 3.0;
  v2 << 1.0, 2.0;

  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v2));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v2));
}

TEST(AgradRevMatrix, test_autodiff_Eigen_RowVectorXd) {
  auto func = [](auto a, auto b) {
    stan::math::check_size_match("test_autodiff_Eigen_RowVectorXd", "a",
                                 a.size(), "b", b.size());
    return (a + b).eval();
  };

  Eigen::RowVectorXd v1(3);
  Eigen::RowVectorXd v2(2);

  v1 << 1.0, 2.0, 3.0;
  v2 << 1.0, 2.0;

  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v2));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v2));
}

TEST(AgradRevMatrix, test_autodiff_Eigen_MatrixXd) {
  auto func = [](auto a, auto b) {
    stan::math::check_matching_dims("test_autodiff_Eigen_MatrixXd", "a", a, "b",
                                    b);
    return (a + b).eval();
  };

  Eigen::MatrixXd v1(3, 2);
  Eigen::MatrixXd v2(2, 2);

  v1 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  v2 << 1.0, 2.0, -1.0, -2.0;

  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v2));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v2));
}

TEST(AgradRevMatrix, test_autodiff_Eigen_MatrixXd_scalar_output) {
  auto func = [](auto a, auto b) {
    stan::math::check_matching_dims("test_autodiff_Eigen_MatrixXd", "a", a, "b",
                                    b);
    return stan::math::sum(a) + stan::math::sum(b);
  };

  Eigen::MatrixXd v1(3, 2);
  Eigen::MatrixXd v2(2, 2);

  v1 << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
  v2 << 1.0, 2.0, -1.0, -2.0;

  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v1, v2));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v1));
  EXPECT_NO_THROW(stan::math::test::test_autodiff(func, v2, v2));
}

struct bad_scalar_func {
  template <typename T>
  T operator()(const T& t) {
    return t;
  }

  double operator()(const double& t) { return t + 1.0; }
};

TEST(AgradRevMatrix, test_autodiff_scalar_bad) {
  EXPECT_THROW_MSG(stan::math::test::test_autodiff(bad_scalar_func{}, 1.0),
                   std::runtime_error,
                   "does not match the value from the reverse mode call");
}

struct bad_std_vector_func {
  template <typename T>
  T operator()(const T& t) {
    return t;
  }

  std::vector<double> operator()(const std::vector<double>& t) {
    std::vector<double> o(t);
    o[0] = o[0] + 1.0;
    return o;
  }
};

TEST(AgradRevMatrix, test_autodiff_std_vector_bad) {
  std::vector<double> v1 = {{1.0, 2.0}};
  EXPECT_THROW_MSG(stan::math::test::test_autodiff(bad_std_vector_func{}, v1),
                   std::runtime_error,
                   "does not match the value from the reverse mode call");
}

struct bad_eigen_vector_func {
  template <typename T>
  T operator()(const T& t) {
    return t;
  }

  Eigen::VectorXd operator()(const Eigen::VectorXd& t) { return t - t; }
};

TEST(AgradRevMatrix, test_autodiff_eigen_vector_bad) {
  Eigen::VectorXd v1(3);

  v1 << 1.0, 2.0, 3.0;

  EXPECT_THROW_MSG(stan::math::test::test_autodiff(bad_eigen_vector_func{}, v1),
                   std::runtime_error,
                   "does not match the value from the reverse mode call");
}

struct bad_eigen_row_vector_func {
  template <typename T>
  T operator()(const T& t) {
    return t;
  }

  Eigen::RowVectorXd operator()(const Eigen::RowVectorXd& t) { return t - t; }
};

TEST(AgradRevMatrix, test_autodiff_eigen_row_vector_bad) {
  Eigen::RowVectorXd v1(3);

  v1 << 1.0, 2.0, 3.0;

  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(bad_eigen_row_vector_func{}, v1),
      std::runtime_error,
      "does not match the value from the reverse mode call");
}

struct bad_eigen_matrix_func {
  template <typename T>
  T operator()(const T& t) {
    return t;
  }

  Eigen::MatrixXd operator()(const Eigen::MatrixXd& t) { return t - t; }
};

TEST(AgradRevMatrix, test_autodiff_eigen_matrix_bad) {
  Eigen::MatrixXd m1(2, 2);

  m1 << 1.0, 2.0, 3.0, 4.0;

  EXPECT_THROW_MSG(stan::math::test::test_autodiff(bad_eigen_matrix_func{}, m1),
                   std::runtime_error,
                   "does not match the value from the reverse mode call");
}

struct bad_autodiff_scalar_func {
  stan::math::var operator()(const stan::math::var& t) {
    return stan::math::to_var(stan::math::value_of(t));
  }

  double operator()(const double& t) { return t; }
};

TEST(AgradRevMatrix, test_autodiff_scalar_bad_autodiff) {
  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(bad_autodiff_scalar_func{}, 1.0),
      std::runtime_error,
      "Jacobian element of the finite difference approximation");
}

struct bad_autodiff_std_vector_func {
  template <typename T>
  T operator()(const T& t) {
    return stan::math::to_var(stan::math::value_of(t));
  }

  std::vector<double> operator()(const std::vector<double>& t) { return t; }
};

TEST(AgradRevMatrix, test_autodiff_std_vector_bad_autodiff) {
  std::vector<double> v1 = {{1.0, 2.0}};
  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(bad_autodiff_std_vector_func{}, v1),
      std::runtime_error,
      "Jacobian element of the finite difference approximation");
}

struct bad_autodiff_eigen_vector_func {
  template <typename T>
  T operator()(const T& t) {
    return stan::math::to_var(stan::math::value_of(t));
  }

  Eigen::VectorXd operator()(const Eigen::VectorXd& t) { return t; }
};

TEST(AgradRevMatrix, test_autodiff_eigen_vector_bad_autodiff) {
  Eigen::VectorXd v1(3);

  v1 << 1.0, 2.0, 3.0;

  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(bad_autodiff_eigen_vector_func{}, v1),
      std::runtime_error,
      "Jacobian element of the finite difference approximation");
}

struct bad_autodiff_eigen_row_vector_func {
  template <typename T>
  T operator()(const T& t) {
    return stan::math::to_var(stan::math::value_of(t));
  }

  Eigen::RowVectorXd operator()(const Eigen::RowVectorXd& t) { return t; }
};

TEST(AgradRevMatrix, test_autodiff_eigen_row_vector_bad_autodiff) {
  Eigen::RowVectorXd v1(3);

  v1 << 1.0, 2.0, 3.0;

  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(bad_autodiff_eigen_row_vector_func{}, v1),
      std::runtime_error,
      "Jacobian element of the finite difference approximation");
}

struct bad_autodiff_eigen_matrix_func {
  template <typename T>
  T operator()(const T& t) {
    return stan::math::to_var(stan::math::value_of(t));
  }

  Eigen::MatrixXd operator()(const Eigen::MatrixXd& t) { return t; }
};

TEST(AgradRevMatrix, test_autodiff_eigen_matrix_bad_autodiff) {
  Eigen::MatrixXd m1(2, 2);

  m1 << 1.0, 2.0, 3.0, 4.0;

  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(bad_autodiff_eigen_matrix_func{}, m1),
      std::runtime_error,
      "Jacobian element of the finite difference approximation");
}

template <typename T, typename E>
struct single_error_helper {
  E e_;

  single_error_helper(E e) : e_(e) {}

  template <typename R>
  R operator()(R a) {
    return a;
  }

  T operator()(T a) { throw e_; }
};

template <typename E1, typename E2>
struct double_error_helper {
  E1 e1_;
  E2 e2_;

  double_error_helper(E1 e1, E2 e2) : e1_(e1), e2_(e2) {}

  template <typename R>
  R operator()(R a) {
    throw e2_;
  }

  double operator()(double a) { throw e1_; }
};

class fake_exception {};

template <typename E>
void test_single_error(E e) {
  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(single_error_helper<double, E>(e), 1.0),
      std::runtime_error, "var version threw nothing");
  EXPECT_THROW_MSG(stan::math::test::test_autodiff(
                       single_error_helper<stan::math::var, E>(e), 1.0),
                   std::runtime_error,
                   "prim version of function threw no exception but");
}

template <typename E>
void test_errors_match(E e) {
  EXPECT_NO_THROW(
      stan::math::test::test_autodiff(double_error_helper<E, E>(e, e), 1.0));
}

TEST(AgradRevMatrix, test_autodiff_single_errors) {
  test_single_error(std::runtime_error(""));
  test_single_error(std::domain_error(""));
  test_single_error(std::logic_error(""));
  test_single_error(std::invalid_argument(""));
  test_single_error(std::out_of_range(""));
  test_single_error(std::system_error(std::error_code()));
  test_single_error(std::exception());
  test_single_error(fake_exception());
}

TEST(AgradRevMatrix, test_autodiff_errors_match) {
  test_errors_match(std::runtime_error(""));
  test_errors_match(std::domain_error(""));
  test_errors_match(std::logic_error(""));
  test_errors_match(std::invalid_argument(""));
  test_errors_match(std::out_of_range(""));
  test_errors_match(std::system_error(std::error_code()));
  test_errors_match(std::exception());
}

template <typename E>
void test_errors_different(E e1, E e2) {}

template <typename E1, typename E2>
void test_errors_different(E1 e1, E2 e2) {
  EXPECT_THROW_MSG(
      stan::math::test::test_autodiff(double_error_helper<E1, E2>(e1, e2), 1.0),
      std::runtime_error, "threw different error");
}

TEST(AgradRevMatrix, text_autodiff_errors_different) {
  stan::math::call_all_argument_combos(
      [](auto e1, auto e2) {
        test_errors_different(e1, e2);
        return 0;
      },
      std::make_tuple(std::domain_error(""), std::invalid_argument(""),
                      std::out_of_range(""),
                      std::system_error(std::error_code())),
      std::make_tuple(std::domain_error(""), std::invalid_argument(""),
                      std::out_of_range(""),
                      std::system_error(std::error_code())));
}
