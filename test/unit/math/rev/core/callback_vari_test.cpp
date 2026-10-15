#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST_F(AgradRev, RevCBack_callback_vari_scalar_test) {
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

TEST_F(AgradRev, RevCBack_callback_vari_const_scalar_compile_test) {
  stan::math::var a = 1;

  const double& a_val = a.val();

  stan::math::var b = stan::math::make_callback_vari(
      a_val, [a](const auto& vi) mutable { a.adj() += vi.adj(); });
  EXPECT_FLOAT_EQ(a.val(), b.val());
}

TEST_F(AgradRev, RevCBack_callback_vari_eigen_test) {
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

TEST_F(AgradRev, RevCBack_callback_vari_over_aligned_capture) {
  // a callback stores its captures inside the vari allocated on the arena
  const Eigen::Vector2d m16(1.0, 2.0);
  const Eigen::Matrix4d m32 = Eigen::Matrix4d::Identity();
  stan::math::var a = 1;
  bool aligned16 = false;
  bool aligned32 = false;
  // leave the next free byte of the arena 8 bytes past a 16 byte boundary
  auto& arena = stan::math::ChainableStack::instance_->memalloc_;
  while (!stan::math::is_aligned(arena.alloc(8), 16U)) {
  }
  stan::math::var b = stan::math::make_callback_vari(
      a.val(), [a, m16, &aligned16](const auto& vi) mutable {
        aligned16
            = stan::math::is_aligned(m16.data(), alignof(Eigen::Vector2d));
        a.adj() += vi.adj_ * m16.sum();
      });
  while (!stan::math::is_aligned(arena.alloc(8), 16U)) {
  }
  stan::math::var c = stan::math::make_callback_vari(
      b.val(), [b, m32, &aligned32](const auto& vi) mutable {
        aligned32
            = stan::math::is_aligned(m32.data(), alignof(Eigen::Matrix4d));
        b.adj() += vi.adj_ * m32.trace();
      });
  c.grad();
  EXPECT_TRUE(aligned16);
  EXPECT_TRUE(aligned32);
  EXPECT_FLOAT_EQ(a.adj(), 12.0);
}

TEST_F(AgradRev, RevCBack_make_callback_var_scalar_test) {
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

TEST_F(AgradRev, RevCBack_make_callback_var_eigen_test) {
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
