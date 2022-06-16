#include <stan/math/rev.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRev, value_of_rec) {
  using stan::math::value_of_rec;
  using stan::math::var;

  var v_a(5.0);

  EXPECT_FLOAT_EQ(5.0, value_of_rec(v_a));
}

TEST(MathMatrixRevArr, value_of_rec) {
  using stan::math::value_of_rec;
  using stan::math::var;
  using std::vector;

  vector<double> a_vals;
  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;
  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  vector<var> a;
  a = stan::math::to_var(a_vals);
  vector<var> b;
  b = stan::math::to_var(b_vals);

  vector<double> d_a = value_of_rec(a);
  vector<double> d_b = value_of_rec(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val(), d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val(), d_a[i]);
}

TEST(AgradMatrixRev, value_of_rec) {
  using stan::math::value_of_rec;
  using stan::math::var;
  using std::vector;

  vector<double> a_vals;

  for (size_t i = 0; i < 10; ++i)
    a_vals.push_back(i + 1);

  vector<double> b_vals;

  for (size_t i = 10; i < 15; ++i)
    b_vals.push_back(i + 1);

  Eigen::Matrix<double, 2, 5> a;
  ::fill(a_vals, a);
  Eigen::Matrix<double, 5, 1> b;
  ::fill(b_vals, b);

  Eigen::Matrix<var, 2, 5> v_a;
  ::fill(a_vals, v_a);
  Eigen::Matrix<var, 5, 1> v_b;
  ::fill(b_vals, v_b);

  Eigen::MatrixXd d_v_a = value_of_rec(v_a);
  Eigen::MatrixXd d_v_b = value_of_rec(v_b);

  for (Eigen::Index i = 0; i < 5; ++i) {
    EXPECT_FLOAT_EQ(b(i), d_v_b(i));
  }

  for (Eigen::Index i = 0; i < 2; ++i)
    for (Eigen::Index j = 0; j < 5; ++j) {
      EXPECT_FLOAT_EQ(a(i, j), d_v_a(i, j));
    }
}

TEST(AgradMatrixRev, value_of_rec_expression) {
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::value_of;
  using stan::math::var;
  Matrix<var, -1, -1> a = MatrixXd::Random(7, 4);
  MatrixXd res = value_of_rec(2 * a);
  MatrixXd correct = 2 * value_of_rec(a);

  EXPECT_MATRIX_NEAR(res, correct, 1e-10);
}

TEST(AgradMatrixRev, value_of_rec_matrix_rvalue) {
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::value_of;
  using stan::math::var;
  Matrix<var, -1, -1> a = MatrixXd::Random(7, 4);
  MatrixXd correct = value_of_rec(a);
  const auto& tmp = value_of_rec(std::move(a));
  // we are expecting an expression, not a plain matrix
  EXPECT_FALSE((std::is_same<std::decay_t<decltype(tmp)>,
                             stan::plain_type_t<decltype(tmp)>>::value));
  a.setZero();
  MatrixXd res = tmp;

  EXPECT_MATRIX_NEAR(res, correct, 1e-10);
}
