#include <test/unit/math/expr_tests.hpp>
#include <gtest/gtest-spi.h>

namespace stan {
namespace test {



template <typename T>
auto fail_single_call_fun(T&& x) {
  Eigen::Matrix<scalar_type_t<T>, -1, 1> res(x.size());
  for (int i = 0; i < x.size(); ++i) {
    res[i] = x[i] * x[i];
  }
  return res;
}
TEST(ExpressionTest, simple_single_fail) {
  auto f_fail = [](auto&& x) {return fail_single_call_fun(x);};
  Eigen::VectorXd x = Eigen::VectorXd::Random(2);
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<double>(f_fail, x), "");
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<stan::math::var>(f_fail, x), "");
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<stan::math::fvar<double>>(f_fail, x), "");
}
template <typename T>
auto pass_single_call_fun(T&& x) {
  auto&& x_ref = stan::math::to_ref(x);
  Eigen::Matrix<scalar_type_t<T>, -1, 1> res(x_ref.size());
  for (int i = 0; i < x_ref.size(); ++i) {
    res[i] = x_ref[i];
  }
  return res;
}
TEST(ExpressionTest, simple_single_pass) {
  auto f_pass = [](auto&& x) {return pass_single_call_fun(x);};
  Eigen::VectorXd x = Eigen::VectorXd::Random(2);
  stan::test::check_expr_test(f_pass, x);
}

template <typename T1, typename T2>
auto fail_double_call_fun(T1&& x1, T2&& x2) {
  Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(x1.size());
  for (int i = 0; i < x1.size(); ++i) {
    res[i] = x1[i] * x2[i] + x1[i] * stan::math::sum(x2);
  }
  return res;
}

TEST(ExpressionTest, simple_double_fail) {
  auto f_fail = [](auto&& x1, auto&& x2) {return fail_double_call_fun(x1, x2);};
  Eigen::VectorXd x1 = Eigen::VectorXd::Random(2);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(2);
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x1, x2), ""), "6 failures");
  std::vector<double> x_vec{1, 2};
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<double>(f_fail, x1, x_vec), "");
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<stan::math::var>(f_fail, x1, x_vec), "");
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<stan::math::fvar<double>>(f_fail, x1, x_vec), "");
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<double>(f_fail, x_vec, x1), "");
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<stan::math::var>(f_fail, x_vec, x1), "");
  EXPECT_NONFATAL_FAILURE(stan::test::internal::check_expr_test<stan::math::fvar<double>>(f_fail, x_vec, x1), "");
}

template <typename T1, typename T2>
auto pass_double_call_fun(T1&& x1, T2&& x2) {
  auto&& x1_ref = stan::math::to_ref(x1);
  auto&& x2_ref = stan::math::to_ref(x2);
  Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(x1_ref.size());
  for (int i = 0; i < x1_ref.size(); ++i) {
    res[i] = x1_ref[i] * x1_ref[i] + x2_ref[i] * stan::math::sum(x2_ref);
  }
  return res;
}

TEST(ExpressionTest, simple_double_pass) {
  auto f_pass = [](auto&& x1, auto&& x2) {return pass_double_call_fun(x1, x2);};
  Eigen::VectorXd x1 = Eigen::VectorXd::Random(2);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(2);
  stan::test::check_expr_test(f_pass, x1, x2);
  std::vector<double> x_vec{1, 2};
  stan::test::check_expr_test(f_pass, x1, x_vec);
  stan::test::check_expr_test(f_pass, x_vec, x1);
}

template <typename T1, typename T2, typename T3>
auto fail_triple_call_fun(T1&& x1, T2&& x2, T3&& x3) {
  Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(x1.size());
  for (int i = 0; i < x1.size(); ++i) {
    res[i] = x1[i] * x1[i] + x2[i] * x2[i] + stan::math::sum(x3);
  }
  return res;
}

TEST(ExpressionTest, simple_triple_fail) {
  auto f_fail = [](auto&& x1, auto&& x2, auto&& x3) {return fail_triple_call_fun(x1, x2, x3);};
  Eigen::VectorXd x1 = Eigen::VectorXd::Random(2);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(2);
  Eigen::VectorXd x3 = Eigen::VectorXd::Random(2);
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x1, x2, x3), ""), "9 failures");
  std::vector<double> x_vec{1, 2};
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x_vec, x2, x3), ""), "6 failures");
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x1, x_vec, x3), ""), "6 failures");
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x1, x2, x_vec), ""), "6 failures");
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x1, x_vec, x_vec), ""), "3 failures");
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x_vec, x2, x_vec), ""), "3 failures");
  EXPECT_NONFATAL_FAILURE(EXPECT_NONFATAL_FAILURE(stan::test::check_expr_test(f_fail, x_vec, x_vec, x3), ""), "3 failures");
}

template <typename T1, typename T2, typename T3>
auto pass_triple_call_fun(T1&& x1, T2&& x2, T3&& x3) {
  auto&& x1_ref = stan::math::to_ref(x1);
  auto&& x2_ref = stan::math::to_ref(x2);
  auto&& x3_ref = stan::math::to_ref(x3);
  Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> res(x1_ref.size());
  for (int i = 0; i < x1_ref.size(); ++i) {
    res[i] = x1_ref[i] * x1_ref[i] + x2_ref[i] * x2_ref[i] + x3_ref[i] * x3_ref[i];
  }
  return res;
}

TEST(ExpressionTest, simple_triple_pass) {
  auto f_pass = [](auto&& x1, auto&& x2, auto&& x3) {return pass_triple_call_fun(x1, x2, x3);};
  Eigen::VectorXd x1 = Eigen::VectorXd::Random(2);
  Eigen::VectorXd x2 = Eigen::VectorXd::Random(2);
  Eigen::VectorXd x3 = Eigen::VectorXd::Random(2);
  stan::test::check_expr_test(f_pass, x1, x2, x3);
  std::vector<double> x_vec{1, 2};
  stan::test::check_expr_test(f_pass, x1, x2, x_vec);
  stan::test::check_expr_test(f_pass, x1, x_vec, x3);
  stan::test::check_expr_test(f_pass, x_vec, x2, x3);
  stan::test::check_expr_test(f_pass, x_vec, x_vec, x3);
}



}
}
