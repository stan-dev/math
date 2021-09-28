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

  EXPECT_SAME_TYPE(var, decltype(av));
  EXPECT_SAME_TYPE(var_mat, decltype(bv));
  EXPECT_SAME_TYPE(var_vec, decltype(cv));
  EXPECT_SAME_TYPE(var_row_vec, decltype(dv));

  EXPECT_SAME_TYPE(var_mat, decltype(ev));
  EXPECT_SAME_TYPE(var_vec, decltype(fv));
  EXPECT_SAME_TYPE(var_row_vec, decltype(gv));

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, to_var_value_vec_types) {
  using stan::math::to_var_value;
  using stan::math::var;
  using stan::math::var_value;

  using mat_var_svec
      = std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>>;
  using vec_var_svec = std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>>;
  using row_vec_var_svec = std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>>;

  using var_mat_svec = std::vector<var_value<Eigen::MatrixXd>>;
  using var_vec_svec = std::vector<var_value<Eigen::VectorXd>>;
  using var_row_vec_svec = std::vector<var_value<Eigen::RowVectorXd>>;

  using var_svec = std::vector<var>;

  var_svec a(1, 2.0);
  mat_var_svec b(1, Eigen::MatrixXd(2, 2));
  vec_var_svec c(1, Eigen::VectorXd(2));
  row_vec_var_svec d(1, Eigen::RowVectorXd(2));

  var_mat_svec e(1, Eigen::MatrixXd(2, 2));
  var_vec_svec f(1, Eigen::VectorXd(2));
  var_row_vec_svec g(1, Eigen::RowVectorXd(2));

  auto av = to_var_value(a);
  auto bv = to_var_value(b);
  auto cv = to_var_value(c);
  auto dv = to_var_value(d);

  auto ev = to_var_value(e);
  auto fv = to_var_value(f);
  auto gv = to_var_value(g);

  EXPECT_SAME_TYPE(var_svec, decltype(av));
  EXPECT_SAME_TYPE(var_mat_svec, decltype(bv));
  EXPECT_SAME_TYPE(var_vec_svec, decltype(cv));
  EXPECT_SAME_TYPE(var_row_vec_svec, decltype(dv));

  EXPECT_SAME_TYPE(var_mat_svec, decltype(ev));
  EXPECT_SAME_TYPE(var_vec_svec, decltype(fv));
  EXPECT_SAME_TYPE(var_row_vec_svec, decltype(gv));

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

TEST(AgradRevMatrix, to_var_value_matrix_svec_test) {
  std::vector<stan::math::promote_scalar_t<stan::math::var, Eigen::MatrixXd>>
      matvars;

  for (size_t i = 0; i < 2; ++i) {
    matvars.push_back(Eigen::MatrixXd::Random(2, 2));
  }

  auto varmats = stan::math::to_var_value(matvars);

  for (size_t i = 0; i < 2; ++i) {
    EXPECT_MATRIX_EQ(matvars[i].val(), varmats[i].val());
    stan::math::set_zero_all_adjoints();
    varmats[i].adj() = Eigen::MatrixXd::Random(2, 2);
    stan::math::grad();
    EXPECT_MATRIX_EQ(matvars[i].adj(), varmats[i].adj());
  }
}

TEST(AgradRevMatrix, to_var_value_vector_svec_test) {
  std::vector<stan::math::promote_scalar_t<stan::math::var, Eigen::VectorXd>>
      matvars;

  for (size_t i = 0; i < 2; ++i) {
    matvars.push_back(Eigen::VectorXd::Random(2));
  }

  auto varmats = stan::math::to_var_value(matvars);

  for (size_t i = 0; i < 2; ++i) {
    EXPECT_MATRIX_EQ(matvars[i].val(), varmats[i].val());
    stan::math::set_zero_all_adjoints();
    varmats[i].adj() = Eigen::VectorXd::Random(2);
    stan::math::grad();
    EXPECT_MATRIX_EQ(matvars[i].adj(), varmats[i].adj());
  }
}

TEST(AgradRevMatrix, to_var_value_row_vector_svec_test) {
  std::vector<stan::math::promote_scalar_t<stan::math::var, Eigen::RowVectorXd>>
      matvars;

  for (size_t i = 0; i < 2; ++i) {
    matvars.push_back(Eigen::RowVectorXd::Random(2));
  }

  auto varmats = stan::math::to_var_value(matvars);

  for (size_t i = 0; i < 2; ++i) {
    EXPECT_MATRIX_EQ(matvars[i].val(), varmats[i].val());
    stan::math::set_zero_all_adjoints();
    varmats[i].adj() = Eigen::RowVectorXd::Random(2);
    stan::math::grad();
    EXPECT_MATRIX_EQ(matvars[i].adj(), varmats[i].adj());
  }
}
