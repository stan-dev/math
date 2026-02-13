#include <stan/math/rev.hpp>
#include <test/unit/math/prim/util.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <utility>
#include <vector>

TEST_F(AgradRev, RevFun_adjoint_of_var_scalar_lvalue) {
  stan::math::var x = 2.0;
  auto&& x_adj = stan::math::adjoint_of(x);
  EXPECT_SAME_TYPE(double&, decltype(x_adj));
  x_adj = 1.75;
  EXPECT_FLOAT_EQ(x.adj(), 1.75);
}

TEST_F(AgradRev, RevFun_adjoint_of_var_scalar_rvalue) {
  auto&& x_adj = stan::math::adjoint_of(stan::math::var(2.0));
  EXPECT_SAME_TYPE(double&, decltype(x_adj));

  x_adj = 3.25;
  EXPECT_FLOAT_EQ(x_adj, 3.25);
}

TEST_F(AgradRev, RevFun_adjoint_of_var_value_matrix_lvalue) {
  using stan::math::var_value;
  var_value<Eigen::MatrixXd> x(Eigen::MatrixXd::Zero(2, 2));
  auto&& x_adj = stan::math::adjoint_of(x);
  EXPECT_SAME_TYPE(decltype(x.adj()), decltype(x_adj));

  Eigen::MatrixXd expected(2, 2);
  expected << 1.0, 2.0, 3.0, 4.0;
  x_adj = expected;
  EXPECT_MATRIX_EQ(x.adj(), expected);
}

TEST_F(AgradRev, RevFun_adjoint_of_matrix_var_lvalue) {
  using matrix_v = Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic>;
  matrix_v x = Eigen::MatrixXd::Zero(2, 3).template cast<stan::math::var>();
  auto&& x_adj = stan::math::adjoint_of(x);

  Eigen::MatrixXd expected(2, 3);
  expected << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;
  x_adj = expected;

  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      EXPECT_FLOAT_EQ(x(i, j).adj(), expected(i, j));
    }
  }
}

TEST_F(AgradRev, RevFun_adjoint_of_accepts_eigen_var_value_expression_rvalue) {
  using stan::math::set_zero_all_adjoints;
  using stan::math::var_value;

  var_value<Eigen::MatrixXd> a(Eigen::MatrixXd::Random(2, 2));
  var_value<Eigen::MatrixXd> b(Eigen::MatrixXd::Random(2, 2));
  auto a_plus_b = a + b;
  decltype(auto) expr_adj = stan::math::adjoint_of(a_plus_b);
  auto sum_var = stan::math::sum(a_plus_b);
  stan::math::grad(sum_var.vi_);
  Eigen::MatrixXd expected = Eigen::MatrixXd::Constant(2, 2, 1);
  EXPECT_MATRIX_EQ(expr_adj, expected);
  EXPECT_MATRIX_EQ(a.adj(), expected);
  EXPECT_MATRIX_EQ(b.adj(), expected);
}

TEST_F(AgradRev, RevFun_adjoint_of_std_vector_var_lvalue) {
  std::vector<stan::math::var> x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  const auto x_size_before = x.size();

  auto&& x_adj = stan::math::adjoint_of(x);
  EXPECT_EQ(x.size(), x_size_before);

  Eigen::VectorXd expected(3);
  expected << 2.0, 4.0, 6.0;
  x_adj = expected;

  EXPECT_FLOAT_EQ(x[0].adj(), 2.0);
  EXPECT_FLOAT_EQ(x[1].adj(), 4.0);
  EXPECT_FLOAT_EQ(x[2].adj(), 6.0);
}

TEST_F(AgradRev, RevFun_adjoint_of_std_vector_var_rvalue) {
  auto x_adj = stan::math::adjoint_of(
      std::vector<stan::math::var>{stan::math::var(1.0), stan::math::var(2.0)});

  EXPECT_EQ(x_adj.size(), 2);
  Eigen::VectorXd expected(2);
  expected << 5.0, 7.0;
  x_adj = expected;
  EXPECT_MATRIX_EQ(x_adj, expected);
}

TEST_F(AgradRev, RevFun_adjoint_of_nested_std_vector_var_lvalue) {
  std::vector<std::vector<stan::math::var>> x{
      {stan::math::var(1.0), stan::math::var(2.0)},
      {stan::math::var(3.0)}};

  const auto outer_size_before = x.size();
  const auto first_inner_size_before = x[0].size();
  const auto second_inner_size_before = x[1].size();
  auto&& x_adj = stan::math::adjoint_of(x);
  EXPECT_EQ(x.size(), outer_size_before);
  EXPECT_EQ(x[0].size(), first_inner_size_before);
  EXPECT_EQ(x[1].size(), second_inner_size_before);

  Eigen::VectorXd expected0(2);
  expected0 << 1.5, 2.5;
  Eigen::VectorXd expected1(1);
  expected1 << 3.5;
  x_adj[0] = expected0;
  x_adj[1] = expected1;

  EXPECT_FLOAT_EQ(x[0][0].adj(), 1.5);
  EXPECT_FLOAT_EQ(x[0][1].adj(), 2.5);
  EXPECT_FLOAT_EQ(x[1][0].adj(), 3.5);
}

TEST_F(AgradRev, RevFun_adjoint_of_nested_std_vector_var_rvalue) {
  auto x_adj = stan::math::adjoint_of(
      std::vector<std::vector<stan::math::var>>{
          {stan::math::var(1.0), stan::math::var(2.0)},
          {stan::math::var(3.0), stan::math::var(4.0), stan::math::var(5.0)}});
  ASSERT_EQ(x_adj.size(), 2);
  EXPECT_EQ(x_adj[0].size(), 2);
  EXPECT_EQ(x_adj[1].size(), 3);
}

TEST_F(AgradRev, RevFun_adjoint_of_tuple_mixed_inputs) {
  using stan::test::is_lvalue_ref_element_v;

  stan::math::var scalar = 1.0;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> mat
      = Eigen::MatrixXd::Zero(2, 2).template cast<stan::math::var>();
  std::vector<stan::math::var> vec{stan::math::var(0.0), stan::math::var(0.0),
                                   stan::math::var(0.0)};

  auto out = stan::math::adjoint_of(std::forward_as_tuple(scalar, mat, vec));
  using out_t = decltype(out);
  EXPECT_EQ(std::tuple_size_v<out_t>, 3);
  EXPECT_TRUE((is_lvalue_ref_element_v<0, out_t>));

  std::get<0>(out) = 7.0;
  Eigen::MatrixXd mat_expected = Eigen::MatrixXd::Constant(2, 2, 3.0);
  std::get<1>(out) = mat_expected;
  Eigen::VectorXd vec_expected = Eigen::VectorXd::Constant(3, 4.0);
  std::get<2>(out) = vec_expected;

  EXPECT_FLOAT_EQ(scalar.adj(), 7.0);
  EXPECT_MATRIX_EQ(mat.adj(), mat_expected);
  EXPECT_FLOAT_EQ(vec[0].adj(), 4.0);
  EXPECT_FLOAT_EQ(vec[1].adj(), 4.0);
  EXPECT_FLOAT_EQ(vec[2].adj(), 4.0);
}

TEST_F(AgradRev, RevFun_adjoint_of_empty_tuple) {
  auto out = stan::math::adjoint_of(std::tuple<>{});
  EXPECT_EQ(std::tuple_size_v<decltype(out)>, 0);
}
