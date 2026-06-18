#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>

// Reproduces a compile error: calling `transpose` on a `var_value` that holds
// a non-plain Eigen expression (e.g. a row/block of a `var_value` matrix). The
// resulting `var_value<Eigen::Transpose<...>>` is not a plain type, so
// `make_holder` tries to wrap the `var_value` in an Eigen `Holder` expression,
// which fails because `var_value` has no `StorageKind`.
TEST_F(AgradRev, RevMatrix_transpose_var_block) {
  using stan::math::transpose;
  using stan::math::var_value;

  Eigen::MatrixXd m(2, 3);
  m << 1, 2, 3, 4, 5, 6;
  var_value<Eigen::MatrixXd> m_vv(m);

  // m_vv.row(0) is a var_value<Eigen::Block<Eigen::Map<...>>>
  auto rt = transpose(m_vv.row(0));

  EXPECT_EQ(3, rt.rows());
  EXPECT_EQ(1, rt.cols());
  EXPECT_FLOAT_EQ(1, rt.val()(0));
  EXPECT_FLOAT_EQ(2, rt.val()(1));
  EXPECT_FLOAT_EQ(3, rt.val()(2));
}
