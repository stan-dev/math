#include <stan/math/fwd.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

namespace {
// Change the value of x, then check that the view of the values shows the
// change. A copy of the values would not show it.
template <typename T>
void expect_view_reads_fvar(T& x) {
  using stan::math::as_value_column_array_or_scalar_ref;
  decltype(auto) x_val = as_value_column_array_or_scalar_ref(x);
  using T_val = std::decay_t<decltype(x_val)>;
  EXPECT_EQ(1, T_val::ColsAtCompileTime);
  EXPECT_FLOAT_EQ(x[1].val_, x_val[1]);
  x[1].val_ = 7;
  EXPECT_FLOAT_EQ(7, x_val[1]);
}
}  // namespace

TEST(MathFunFwd, as_value_column_array_or_scalar_ref_scalar) {
  using stan::math::as_value_column_array_or_scalar_ref;
  stan::math::fvar<double> a(2.5, 1.0);
  double a_val = as_value_column_array_or_scalar_ref(a);
  EXPECT_EQ(2.5, a_val);
}

TEST(MathFunFwd, as_value_column_array_or_scalar_ref_reads_fvar) {
  using stan::math::fvar;
  Eigen::Matrix<fvar<double>, -1, 1> a = Eigen::VectorXd::Ones(3);
  Eigen::Matrix<fvar<double>, 1, -1> b = Eigen::RowVectorXd::Ones(3);
  std::vector<fvar<double>> c(3, fvar<double>(1.0));
  expect_view_reads_fvar(a);
  expect_view_reads_fvar(b);
  expect_view_reads_fvar(c);
}
