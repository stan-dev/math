#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <vector>
#include <tuple>
#include <gtest/gtest.h>

TEST(AgradRevZero, zero_arithmetic) {
  int a = 1.0;
  double b = 2;
  std::vector<int> va(5, a);
  std::vector<double> vb(5, b);
  Eigen::VectorXd c = Eigen::VectorXd::Random(5);
  Eigen::RowVectorXd d = Eigen::RowVectorXd::Random(5);
  Eigen::MatrixXd e = Eigen::MatrixXd::Random(5, 5);
  std::vector<std::vector<int>> vva(5, va);
  std::vector<std::vector<double>> vvb(5, vb);
  std::vector<Eigen::VectorXd> vc(5, c);
  std::vector<Eigen::RowVectorXd> vd(5, d);
  std::vector<Eigen::MatrixXd> ve(5, e);

  stan::math::zero_adjoints(a);
  stan::math::zero_adjoints(b);
  stan::math::zero_adjoints(va);
  stan::math::zero_adjoints(vb);
  stan::math::zero_adjoints(c);
  stan::math::zero_adjoints(d);
  stan::math::zero_adjoints(e);
  stan::math::zero_adjoints(vva);
  stan::math::zero_adjoints(vvb);
  stan::math::zero_adjoints(vc);
  stan::math::zero_adjoints(vd);
  stan::math::zero_adjoints(ve);
  stan::math::for_each(
      [](auto&& x) { stan::math::zero_adjoints(x); },
      std::forward_as_tuple(a, b, va, vb, c, d, e, vva, vvb, vc, vd, ve));
}

TEST(AgradRevZero, zero_var) {
  using stan::math::var;
  using stan::math::vari;

  var a(5.0);
  a.vi_->adj_ = 2.0;

  stan::math::zero_adjoints(a);
  EXPECT_FLOAT_EQ(a.vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_std_vector_var) {
  using stan::math::var;
  using stan::math::vari;

  std::vector<var> a(5, 0.0);
  for (size_t i = 0; i < a.size(); ++i)
    a[i].vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(a);
  for (size_t i = 0; i < a.size(); ++i)
    EXPECT_FLOAT_EQ(a[i].vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_vector_var) {
  using stan::math::var;
  using stan::math::vari;

  Eigen::Matrix<var, Eigen::Dynamic, 1> a
      = Eigen::VectorXd::Zero(5).template cast<var>();
  for (size_t i = 0; i < a.size(); ++i)
    a(i).vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(a);
  for (size_t i = 0; i < a.size(); ++i)
    EXPECT_FLOAT_EQ(a(i).vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_row_vector_var) {
  using stan::math::var;
  using stan::math::vari;

  Eigen::Matrix<var, 1, Eigen::Dynamic> a
      = Eigen::RowVectorXd::Zero(5).template cast<var>();
  for (size_t i = 0; i < a.size(); ++i)
    a(i).vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(a);
  for (size_t i = 0; i < a.size(); ++i)
    EXPECT_FLOAT_EQ(a(i).vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_matrix_var) {
  using stan::math::var;
  using stan::math::vari;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a
      = Eigen::MatrixXd::Zero(5, 5).template cast<var>();
  for (size_t i = 0; i < a.size(); ++i)
    a(i).vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(a);
  for (size_t i = 0; i < a.size(); ++i)
    EXPECT_FLOAT_EQ(a(i).vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_std_vector_std_vector_var) {
  using stan::math::var;
  using stan::math::vari;

  std::vector<var> a(5, 0.0);
  std::vector<var> b(5, 1.0);
  std::vector<var> c(5, 2.0);
  std::vector<std::vector<var>> va = {a, b, c};
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      va[i][j].vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(va);
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      EXPECT_FLOAT_EQ(va[i][j].vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_std_vector_vector_var) {
  using stan::math::var;
  using stan::math::vari;

  Eigen::Matrix<var, Eigen::Dynamic, 1> a
      = Eigen::VectorXd::Zero(3).template cast<var>();
  Eigen::Matrix<var, Eigen::Dynamic, 1> b
      = Eigen::VectorXd::Zero(4).template cast<var>();
  Eigen::Matrix<var, Eigen::Dynamic, 1> c
      = Eigen::VectorXd::Zero(5).template cast<var>();
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> va = {a, b, c};
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      va[i](j).vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(va);
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      EXPECT_FLOAT_EQ(va[i](j).vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_std_vector_row_vector_var) {
  using stan::math::var;
  using stan::math::vari;

  Eigen::Matrix<var, 1, Eigen::Dynamic> a
      = Eigen::RowVectorXd::Zero(3).template cast<var>();
  Eigen::Matrix<var, 1, Eigen::Dynamic> b
      = Eigen::RowVectorXd::Zero(4).template cast<var>();
  Eigen::Matrix<var, 1, Eigen::Dynamic> c
      = Eigen::RowVectorXd::Zero(5).template cast<var>();
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> va = {a, b, c};
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      va[i](j).vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(va);
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      EXPECT_FLOAT_EQ(va[i](j).vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_std_vector_matrix_var) {
  using stan::math::var;
  using stan::math::vari;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a
      = Eigen::MatrixXd::Zero(3, 4).template cast<var>();
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> b
      = Eigen::MatrixXd::Zero(4, 5).template cast<var>();
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> c
      = Eigen::MatrixXd::Zero(5, 6).template cast<var>();
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> va
      = {a, b, c};
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      va[i](j).vi_->adj_ = i + 1.0;

  stan::math::zero_adjoints(va);
  for (size_t i = 0; i < va.size(); ++i)
    for (size_t j = 0; j < va[i].size(); ++j)
      EXPECT_FLOAT_EQ(va[i](j).vi_->adj_, 0.0);

  stan::math::recover_memory();
}

TEST(AgradRevZero, zero_multi) {
  using stan::math::var;
  using stan::math::vari;

  int a = 2;
  double b = 3;
  var c(5.0);
  c.vi_->adj_ = 2.0;
  std::vector<var> d(5, 1.0);
  for (size_t i = 0; i < d.size(); ++i)
    d[i].vi_->adj_ = i + 1.0;
  std::vector<int> e(5, 1);
  std::vector<double> f(5, 1.0);

  stan::math::for_each([](auto&& x) { stan::math::zero_adjoints(x); },
                       std::forward_as_tuple(a, b, c, d, e, f));
  EXPECT_FLOAT_EQ(c.vi_->adj_, 0.0);
  for (size_t i = 0; i < d.size(); ++i)
    EXPECT_FLOAT_EQ(d[i].vi_->adj_, 0.0);

  stan::math::recover_memory();
}
