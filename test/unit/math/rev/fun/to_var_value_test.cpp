#include <stan/math/rev.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrix, to_var_value_types) {
  using stan::math::to_var_value;
  using stan::math::var;
  using stan::math::var_value;

  using mat_var = Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>;
  using vec_var = Eigen::Matrix<var, Eigen::Dynamic, 1>;
  using row_vec_var = Eigen::Matrix<var, 1, Eigen::Dynamic>;

  using var_mat = var_value<Eigen::MatrixXd>;
  using var_vec = var_value<Eigen::VectorXd>;
  using var_row_vec = var_value<Eigen::RowVectorXd>;

  var a = 2.0;
  mat_var b = Eigen::MatrixXd(2, 2);
  vec_var c = Eigen::VectorXd(2);
  row_vec_var d = Eigen::RowVectorXd(2);

  var_mat e = Eigen::MatrixXd(2, 2);
  var_vec f = Eigen::VectorXd(2);
  var_row_vec g = Eigen::RowVectorXd(2);

  auto av = to_var_value(a);
  auto bv = to_var_value(b);
  auto cv = to_var_value(c);
  auto dv = to_var_value(d);

  auto ev = to_var_value(e);
  auto fv = to_var_value(f);
  auto gv = to_var_value(g);

  test::expect_same_type<var, decltype(av)>();
  test::expect_same_type<var_mat, decltype(bv)>();
  test::expect_same_type<var_vec, decltype(cv)>();
  test::expect_same_type<var_row_vec, decltype(dv)>();

  test::expect_same_type<var_mat, decltype(ev)>();
  test::expect_same_type<var_vec, decltype(fv)>();
  test::expect_same_type<var_row_vec, decltype(gv)>();

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, to_var_value_if_types) {
  using stan::math::to_var_value_if;
  using stan::math::var;
  using stan::math::var_value;

  using mat_var = Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>;
  using vec_var = Eigen::Matrix<var, Eigen::Dynamic, 1>;
  using row_vec_var = Eigen::Matrix<var, 1, Eigen::Dynamic>;

  using var_mat = var_value<Eigen::MatrixXd>;
  using var_vec = var_value<Eigen::VectorXd>;
  using var_row_vec = var_value<Eigen::RowVectorXd>;

  var a = 2.0;
  mat_var b = Eigen::MatrixXd(2, 2);
  vec_var c = Eigen::VectorXd(2);
  row_vec_var d = Eigen::RowVectorXd(2);

  var_mat e = Eigen::MatrixXd(2, 2);
  var_vec f = Eigen::VectorXd(2);
  var_row_vec g = Eigen::RowVectorXd(2);

  double h = 4.0;
  Eigen::MatrixXd i(2, 2);
  Eigen::VectorXd j(2);
  Eigen::RowVectorXd k(2);

  {
    auto av = to_var_value_if<true>(a);
    auto bv = to_var_value_if<true>(b);
    auto cv = to_var_value_if<true>(c);
    auto dv = to_var_value_if<true>(d);

    auto ev = to_var_value_if<true>(e);
    auto fv = to_var_value_if<true>(f);
    auto gv = to_var_value_if<true>(g);

    test::expect_same_type<var, decltype(av)>();
    test::expect_same_type<var_mat, decltype(bv)>();
    test::expect_same_type<var_vec, decltype(cv)>();
    test::expect_same_type<var_row_vec, decltype(dv)>();

    test::expect_same_type<var_mat, decltype(ev)>();
    test::expect_same_type<var_vec, decltype(fv)>();
    test::expect_same_type<var_row_vec, decltype(gv)>();
  }

  {
    auto av = to_var_value_if<false>(a);
    auto bv = to_var_value_if<false>(b);
    auto cv = to_var_value_if<false>(c);
    auto dv = to_var_value_if<false>(d);

    auto ev = to_var_value_if<false>(e);
    auto fv = to_var_value_if<false>(f);
    auto gv = to_var_value_if<false>(g);

    auto hv = to_var_value_if<false>(h);
    auto iv = to_var_value_if<false>(i);
    auto jv = to_var_value_if<false>(j);

    test::expect_same_type<decltype(a), decltype(av)>();
    test::expect_same_type<decltype(b), decltype(bv)>();
    test::expect_same_type<decltype(c), decltype(cv)>();
    test::expect_same_type<decltype(d), decltype(dv)>();

    test::expect_same_type<decltype(e), decltype(ev)>();
    test::expect_same_type<decltype(f), decltype(fv)>();
    test::expect_same_type<decltype(g), decltype(gv)>();

    test::expect_same_type<decltype(h), decltype(hv)>();
    test::expect_same_type<decltype(i), decltype(iv)>();
    test::expect_same_type<decltype(j), decltype(jv)>();
  }

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, to_var_value_matrix_test) {
  Eigen::MatrixXd val(2, 3);
  val << 1, 2, 3, 4, 5, 6;
  Eigen::MatrixXd adj(2, 3);
  val << 4, 5, 6, 7, 8, 9;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> mat_var = val;
  stan::math::var_value<Eigen::MatrixXd> var_value
      = stan::math::to_var_value(mat_var);
  EXPECT_MATRIX_EQ(var_value.val(), val);
  var_value.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(mat_var.adj(), adj);
}

TEST(AgradRevMatrix, to_var_value_vector_test) {
  Eigen::VectorXd val(3);
  val << 1, 2, 3;
  Eigen::VectorXd adj(3);
  val << 7, 8, 9;
  Eigen::Matrix<stan::math::var, Eigen::Dynamic, 1> mat_var = val;
  stan::math::var_value<Eigen::VectorXd> var_value
      = stan::math::to_var_value(mat_var);
  EXPECT_MATRIX_EQ(var_value.val(), val);
  var_value.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(mat_var.adj(), adj);
}

TEST(AgradRevMatrix, to_var_value_row_vector_test) {
  Eigen::RowVectorXd val(3);
  val << 1, 2, 3;
  Eigen::RowVectorXd adj(3);
  val << 7, 8, 9;
  Eigen::Matrix<stan::math::var, 1, Eigen::Dynamic> mat_var = val;
  stan::math::var_value<Eigen::RowVectorXd> var_value
      = stan::math::to_var_value(mat_var);
  EXPECT_MATRIX_EQ(var_value.val(), val);
  var_value.adj() = adj;
  stan::math::grad();
  EXPECT_MATRIX_EQ(mat_var.adj(), adj);
}
