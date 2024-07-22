#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(AgradRevCBack, callback_vari_scalar_test) {
  stan::math::var a = 1;
  stan::math::var b = 1;

  stan::math::var c = stan::math::make_callback_vari(
      a.val() + b.val() + 3,
      [avi = a.vi_, bvi = b.vi_](const auto& vi) mutable {
        avi->adj_ += vi.adj_;
        bvi->adj_ += vi.adj_ + 1;
      });
  EXPECT_EQ(c.val(), 5);
  EXPECT_EQ(a.adj(), 0);
  EXPECT_EQ(b.adj(), 0);
  c.grad();
  EXPECT_EQ(a.adj(), 1);
  EXPECT_EQ(b.adj(), 2);
}

TEST(AgradRevCBack, callback_vari_const_scalar_compile_test) {
  stan::math::var a = 1;

  const double& a_val = a.val();

  stan::math::var b = stan::math::make_callback_vari(
      a_val, [a](const auto& vi) mutable { a.adj() += vi.adj(); });
  EXPECT_FLOAT_EQ(a.val(), b.val());
}

TEST(AgradRevCBack, callback_vari_eigen_test) {
  Eigen::MatrixXd val(2, 3);
  val << 1, 2, 3, 4, 5, 6;
  stan::math::var_value<Eigen::MatrixXd> a = val;

  stan::math::var_value<Eigen::MatrixXd> b = stan::math::make_callback_vari(
      (a.val().array() + 1).matrix(),
      [avi = a.vi_](const auto& vi) mutable { avi->adj_ += vi.adj_ * 2; });
  EXPECT_MATRIX_EQ(b.val(), (val.array() + 1).matrix());
  EXPECT_MATRIX_EQ(a.adj(), Eigen::MatrixXd::Zero(2, 3));
  stan::math::sum(b).grad();
  EXPECT_MATRIX_EQ(a.adj(), Eigen::MatrixXd::Constant(2, 3, 2));
}

TEST(AgradRevCBack, make_callback_var_scalar_test) {
  stan::math::var a = 1;
  stan::math::var b = 1;

  auto c = stan::math::make_callback_var(
      a.val() + b.val() + 3,
      [avi = a.vi_, bvi = b.vi_](const auto& vi) mutable {
        avi->adj_ += vi.adj_;
        bvi->adj_ += vi.adj_ + 1;
      });
  EXPECT_TRUE((std::is_same<stan::math::var, decltype(c)>::value));
  EXPECT_EQ(c.val(), 5);
  EXPECT_EQ(a.adj(), 0);
  EXPECT_EQ(b.adj(), 0);
  c.grad();
  EXPECT_EQ(a.adj(), 1);
  EXPECT_EQ(b.adj(), 2);
}

TEST(AgradRevCBack, make_callback_var_eigen_test) {
  Eigen::MatrixXd val(2, 3);
  val << 1, 2, 3, 4, 5, 6;
  stan::math::var_value<Eigen::MatrixXd> a = val;

  auto b = stan::math::make_callback_var(
      (a.val().array() + 1).matrix(),
      [avi = a.vi_](const auto& vi) mutable { avi->adj_ += vi.adj_ * 2; });

  EXPECT_TRUE((std::is_same<stan::math::var_value<Eigen::MatrixXd>,
                            decltype(b)>::value));
  EXPECT_MATRIX_EQ(b.val(), (val.array() + 1).matrix());
  EXPECT_MATRIX_EQ(a.adj(), Eigen::MatrixXd::Zero(2, 3));
  stan::math::sum(b).grad();
  EXPECT_MATRIX_EQ(a.adj(), Eigen::MatrixXd::Constant(2, 3, 2));
}
