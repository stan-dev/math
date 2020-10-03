#include <stan/math/rev.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(AgradRev, quad_form_sym_return_types) {
  using stan::math::quad_form_sym;
  using stan::math::var;

  var a = 5.0;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> A(2, 2);
  A << 1.0, 1.0, 1.0, 1.0;

  Eigen::Matrix<var, Eigen::Dynamic, 1> b(2);
  b << 1.0, 1.0;

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> B(2, 1);
  B << 1.0, 1.0;

  EXPECT_TRUE((std::is_same<var, decltype(quad_form_sym(A, b))>::value));
  EXPECT_TRUE((std::is_same<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>,
                            decltype(quad_form_sym(A, B))>::value));
}
