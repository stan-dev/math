#include <stan/math/fwd.hpp>
#include <stan/math/rev.hpp>
#include <test/unit/math/prim/util.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <tuple>
#include <utility>
#include <vector>

TEST_F(AgradRev, FwdFun_tangent_of_fvar_var_scalar_lvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;

  fvar_var x(stan::math::var(2.0), stan::math::var(0.0));
  auto&& x_tan = stan::math::tangent_of(x);
  EXPECT_SAME_TYPE(stan::math::var&, decltype(x_tan));
  x_tan = 1.75;
  EXPECT_FLOAT_EQ(x.d_.val(), 1.75);
}

TEST_F(AgradRev, FwdFun_tangent_of_fvar_var_scalar_rvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;

  auto&& x_tan
      = stan::math::tangent_of(fvar_var(stan::math::var(2.0), stan::math::var(0.0)));
  EXPECT_SAME_TYPE(stan::math::var&, decltype(x_tan));

  x_tan = 3.25;
  EXPECT_FLOAT_EQ(x_tan.val(), 3.25);
}

TEST_F(AgradRev, FwdFun_tangent_of_matrix_fvar_var_lvalue_type_and_values) {
  using fvar_var = stan::math::fvar<stan::math::var>;
  using matrix_fvar_var
      = Eigen::Matrix<fvar_var, Eigen::Dynamic, Eigen::Dynamic>;

  matrix_fvar_var x(2, 2);
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      x(i, j) = fvar_var(stan::math::var(0.0), stan::math::var(0.0));
    }
  }
  decltype(auto) x_tan = stan::math::tangent_of(x);
  EXPECT_SAME_TYPE(decltype(x.d()), decltype(x_tan));

  Eigen::MatrixXd expected(2, 2);
  expected << 1.0, 2.0, 3.0, 4.0;
  x_tan = expected;
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      EXPECT_FLOAT_EQ(x(i, j).d_.val(), expected(i, j));
    }
  }
}

TEST_F(AgradRev, FwdFun_tangent_of_matrix_fvar_var_lvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;
  using matrix_fvar_var
      = Eigen::Matrix<fvar_var, Eigen::Dynamic, Eigen::Dynamic>;

  matrix_fvar_var x(2, 3);
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      x(i, j) = fvar_var(stan::math::var(0.0), stan::math::var(0.0));
    }
  }
  auto&& x_tan = stan::math::tangent_of(x);

  Eigen::MatrixXd expected(2, 3);
  expected << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;
  x_tan = expected;

  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    for (Eigen::Index j = 0; j < x.cols(); ++j) {
      EXPECT_FLOAT_EQ(x(i, j).d_.val(), expected(i, j));
    }
  }
}

TEST_F(AgradRev, FwdFun_tangent_of_accepts_eigen_matrix_fvar_var_expression_rvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;
  using matrix_fvar_var
      = Eigen::Matrix<fvar_var, Eigen::Dynamic, Eigen::Dynamic>;

  matrix_fvar_var a(2, 2);
  matrix_fvar_var b(2, 2);
  for (Eigen::Index i = 0; i < a.rows(); ++i) {
    for (Eigen::Index j = 0; j < a.cols(); ++j) {
      const double aval = static_cast<double>(i + j + 1);
      const double bval = static_cast<double>(2 * i + j + 0.5);
      a(i, j) = fvar_var(stan::math::var(0.0), stan::math::var(aval));
      b(i, j) = fvar_var(stan::math::var(0.0), stan::math::var(bval));
    }
  }
  auto a_plus_b = a + b;
  decltype(auto) expr_tan = stan::math::tangent_of(a_plus_b);
  Eigen::MatrixXd expected(2, 2);
  expected << 1.5, 3.5, 4.5, 6.5;
  EXPECT_MATRIX_EQ(stan::math::value_of(expr_tan), expected);
}

TEST_F(AgradRev, FwdFun_tangent_of_std_vector_fvar_var_lvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;

  std::vector<fvar_var> x(3);
  x[0] = fvar_var(stan::math::var(1.0), stan::math::var(0.0));
  x[1] = fvar_var(stan::math::var(2.0), stan::math::var(0.0));
  x[2] = fvar_var(stan::math::var(3.0), stan::math::var(0.0));
  const auto x_size_before = x.size();

  auto&& x_tan = stan::math::tangent_of(x);
  EXPECT_EQ(x_tan.size(), x_size_before);
  x_tan[0] = 2.0;
  x_tan[1] = 4.0;
  x_tan[2] = 6.0;

  EXPECT_FLOAT_EQ(x[0].d_.val(), 2.0);
  EXPECT_FLOAT_EQ(x[1].d_.val(), 4.0);
  EXPECT_FLOAT_EQ(x[2].d_.val(), 6.0);
}

TEST_F(AgradRev, FwdFun_tangent_of_std_vector_fvar_var_rvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;

  auto x_tan
      = stan::math::tangent_of(std::vector<fvar_var>{
          fvar_var(stan::math::var(1.0), stan::math::var(0.0)),
          fvar_var(stan::math::var(2.0), stan::math::var(0.0))});

  EXPECT_EQ(x_tan.size(), 2);
  x_tan[0] = 5.0;
  x_tan[1] = 7.0;
  EXPECT_FLOAT_EQ(x_tan[0].val(), 5.0);
  EXPECT_FLOAT_EQ(x_tan[1].val(), 7.0);
}

TEST_F(AgradRev, FwdFun_tangent_of_nested_std_vector_fvar_var_lvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;

  std::vector<std::vector<fvar_var>> x{
      {fvar_var(stan::math::var(1.0), stan::math::var(0.0)),
       fvar_var(stan::math::var(2.0), stan::math::var(0.0))},
      {fvar_var(stan::math::var(3.0), stan::math::var(0.0))}};

  const auto outer_size_before = x.size();
  const auto first_inner_size_before = x[0].size();
  const auto second_inner_size_before = x[1].size();
  auto&& x_tan = stan::math::tangent_of(x);
  EXPECT_EQ(x_tan.size(), outer_size_before);
  EXPECT_EQ(x_tan[0].size(), first_inner_size_before);
  EXPECT_EQ(x_tan[1].size(), second_inner_size_before);

  x_tan[0][0] = 1.5;
  x_tan[0][1] = 2.5;
  x_tan[1][0] = 3.5;

  EXPECT_FLOAT_EQ(x[0][0].d_.val(), 1.5);
  EXPECT_FLOAT_EQ(x[0][1].d_.val(), 2.5);
  EXPECT_FLOAT_EQ(x[1][0].d_.val(), 3.5);
}

TEST_F(AgradRev, FwdFun_tangent_of_nested_std_vector_fvar_var_rvalue) {
  using fvar_var = stan::math::fvar<stan::math::var>;

  auto x_tan = stan::math::tangent_of(std::vector<std::vector<fvar_var>>{
      {fvar_var(stan::math::var(1.0), stan::math::var(0.0)),
       fvar_var(stan::math::var(2.0), stan::math::var(0.0))},
      {fvar_var(stan::math::var(3.0), stan::math::var(0.0)),
       fvar_var(stan::math::var(4.0), stan::math::var(0.0)),
       fvar_var(stan::math::var(5.0), stan::math::var(0.0))}});
  ASSERT_EQ(x_tan.size(), 2);
  EXPECT_EQ(x_tan[0].size(), 2);
  EXPECT_EQ(x_tan[1].size(), 3);
}

TEST_F(AgradRev, FwdFun_tangent_of_tuple_mixed_inputs) {
  using fvar_var = stan::math::fvar<stan::math::var>;
  using stan::test::is_lvalue_ref_element_v;

  fvar_var scalar(stan::math::var(1.0), stan::math::var(0.0));
  Eigen::Matrix<fvar_var, Eigen::Dynamic, Eigen::Dynamic> mat(2, 2);
  for (Eigen::Index i = 0; i < mat.rows(); ++i) {
    for (Eigen::Index j = 0; j < mat.cols(); ++j) {
      mat(i, j) = fvar_var(stan::math::var(0.0), stan::math::var(0.0));
    }
  }
  std::vector<fvar_var> vec{
      fvar_var(stan::math::var(0.0), stan::math::var(0.0)),
      fvar_var(stan::math::var(0.0), stan::math::var(0.0)),
      fvar_var(stan::math::var(0.0), stan::math::var(0.0))};

  auto out = stan::math::tangent_of(std::forward_as_tuple(scalar, mat, vec));
  using out_t = decltype(out);
  EXPECT_EQ(std::tuple_size_v<out_t>, 3);
  EXPECT_TRUE((is_lvalue_ref_element_v<0, out_t>));

  std::get<0>(out) = 7.0;
  Eigen::MatrixXd mat_expected = Eigen::MatrixXd::Constant(2, 2, 3.0);
  std::get<1>(out) = mat_expected;
  for (int i = 0; i < 3; ++i) {
    std::get<2>(out)[i] = 4.0;
  }

  EXPECT_FLOAT_EQ(scalar.d_.val(), 7.0);
  for (Eigen::Index i = 0; i < mat.rows(); ++i) {
    for (Eigen::Index j = 0; j < mat.cols(); ++j) {
      EXPECT_FLOAT_EQ(mat(i, j).d_.val(), mat_expected(i, j));
    }
  }
  EXPECT_FLOAT_EQ(vec[0].d_.val(), 4.0);
  EXPECT_FLOAT_EQ(vec[1].d_.val(), 4.0);
  EXPECT_FLOAT_EQ(vec[2].d_.val(), 4.0);
}

TEST_F(AgradRev, FwdFun_tangent_of_empty_tuple) {
  auto out = stan::math::tangent_of(std::tuple<>{});
  EXPECT_EQ(std::tuple_size_v<decltype(out)>, 0);
}
