#include <stan/math/rev.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRev, value_of) {
  using stan::math::value_of;
  using stan::math::var;

  var a = 5.0;
  EXPECT_FLOAT_EQ(5.0, value_of(a));
  // make sure all work together
  EXPECT_FLOAT_EQ(5.0, value_of(5.0));
  EXPECT_FLOAT_EQ(5.0, value_of(5));
}

TEST(MathMatrixRevArr, value_of) {
  using stan::math::value_of;
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

  vector<double> d_a = value_of(a);
  vector<double> d_b = value_of(b);

  for (int i = 0; i < 5; ++i)
    EXPECT_FLOAT_EQ(b[i].val(), d_b[i]);

  for (int i = 0; i < 10; ++i)
    EXPECT_FLOAT_EQ(a[i].val(), d_a[i]);
}

TEST(AgradMatrix, value_of) {
  using stan::math::value_of;
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

  Eigen::MatrixXd d_a = value_of(a);
  Eigen::VectorXd d_b = value_of(b);
  Eigen::MatrixXd d_v_a = value_of(v_a);
  Eigen::MatrixXd d_v_b = value_of(v_b);

  for (Eigen::Index i = 0; i < 5; ++i) {
    EXPECT_FLOAT_EQ(b(i), d_b(i));
    EXPECT_FLOAT_EQ(b(i), d_v_b(i));
  }

  for (Eigen::Index i = 0; i < 2; ++i)
    for (Eigen::Index j = 0; j < 5; ++j) {
      EXPECT_FLOAT_EQ(a(i, j), d_a(i, j));
      EXPECT_FLOAT_EQ(a(i, j), d_v_a(i, j));
    }
}

TEST(AgradMatrix, value_of_vector_of_vectors) {
  using stan::math::var;
  std::vector<var> a(5, 0);
  const std::vector<var> b(5, 0);
  std::vector<std::vector<var>> va(5, a);
  const std::vector<std::vector<var>> vb(5, b);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(va)),
                            std::vector<std::vector<double>>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(vb)),
                            std::vector<std::vector<double>>>::value));

  auto vva = stan::math::value_of(va);
  auto vvb = stan::math::value_of(va);

  for (size_t i = 0; i < va.size(); ++i) {
    for (size_t j = 0; j < va[i].size(); ++j) {
      EXPECT_FLOAT_EQ(vva[i][j], a[j].val());
    }
  }

  for (size_t i = 0; i < vb.size(); ++i) {
    for (size_t j = 0; j < vb[i].size(); ++j) {
      EXPECT_FLOAT_EQ(vvb[i][j], b[j].val());
    }
  }
}

TEST(AgradMatrix, value_of_vector_of_eigen) {
  using stan::math::var;
  Eigen::Matrix<var, Eigen::Dynamic, 1> a
      = Eigen::VectorXd::Random(5).template cast<var>();
  Eigen::Matrix<var, 1, Eigen::Dynamic> b
      = Eigen::RowVectorXd::Random(5).template cast<var>();
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> c
      = Eigen::MatrixXd::Random(5, 5).template cast<var>();
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> va(5, a);
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> vb(5, b);
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> vc(5, c);
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(va)),
                            std::vector<Eigen::VectorXd>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(vb)),
                            std::vector<Eigen::RowVectorXd>>::value));
  EXPECT_TRUE((std::is_same<decltype(stan::math::value_of(vc)),
                            std::vector<Eigen::MatrixXd>>::value));

  auto vva = stan::math::value_of(va);
  auto vvb = stan::math::value_of(vb);
  auto vvc = stan::math::value_of(vc);

  for (size_t i = 0; i < vva.size(); ++i)
    for (size_t j = 0; j < vva[i].size(); ++j)
      EXPECT_FLOAT_EQ(vva[i](j), a(j).val());

  for (size_t i = 0; i < vvb.size(); ++i)
    for (size_t j = 0; j < vva[i].size(); ++j)
      EXPECT_FLOAT_EQ(vvb[i](j), b(j).val());

  for (size_t i = 0; i < vvc.size(); ++i)
    for (size_t j = 0; j < vva[i].size(); ++j)
      EXPECT_FLOAT_EQ(vvc[i](j), c(j).val());
}

TEST(AgradMatrix, value_of_expression) {
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::value_of;
  using stan::math::var;
  Matrix<var, -1, -1> a = MatrixXd::Random(7, 4);
  MatrixXd res = value_of(2 * a);
  MatrixXd correct = 2 * value_of(a);

  EXPECT_MATRIX_NEAR(res, correct, 1e-10);
}

TEST(AgradMatrixRev, value_of_matrix_rvalue) {
  using Eigen::Matrix;
  using Eigen::MatrixXd;
  using stan::math::value_of;
  using stan::math::var;
  Matrix<var, -1, -1> a = MatrixXd::Random(7, 4);
  MatrixXd correct = value_of(a);
  const auto& tmp = value_of(std::move(a));
  // we are expecting an expression, not a plain matrix
  EXPECT_FALSE((std::is_same<std::decay_t<decltype(tmp)>,
                             stan::plain_type_t<decltype(tmp)>>::value));
  a.setZero();
  MatrixXd res = tmp;

  EXPECT_MATRIX_NEAR(res, correct, 1e-10);
}
