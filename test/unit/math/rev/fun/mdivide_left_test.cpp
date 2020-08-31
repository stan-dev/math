#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(AgradRev, mdivide_left_rev) {
  using stan::math::var;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> A(2, 2);
  A << 1.5, 0.7,
    0.3, 1.2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> B(2);
  B << 1.1,
    -1.2;
  
  stan::math::set_zero_all_adjoints();
  var y = sum(mdivide_left(A, B));
  y.grad();

  auto adj = A.adj().eval();

  stan::math::set_zero_all_adjoints();
  var y_ref = sum(multiply(inverse(A), B));
  y_ref.grad();

  EXPECT_FLOAT_EQ(y.val(), y_ref.val());
  EXPECT_FLOAT_EQ(adj(0, 0), A.adj()(0, 0));
}
