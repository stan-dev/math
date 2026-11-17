#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <type_traits>
#include <vector>

namespace {
// Change the value of x through its vari, then check that the view of the
// values shows the change. A copy of the values would not show it.
template <typename T>
void expect_view_reads_vari(T& x) {
  using stan::math::as_value_column_array_or_scalar_ref;
  decltype(auto) x_val = as_value_column_array_or_scalar_ref(x);
  using T_val = std::decay_t<decltype(x_val)>;
  EXPECT_EQ(1, T_val::ColsAtCompileTime);
  EXPECT_FLOAT_EQ(x[1].val(), x_val[1]);
  const_cast<double&>(x[1].vi_->val_) = 7;
  EXPECT_FLOAT_EQ(7, x_val[1]);
}
}  // namespace

TEST_F(AgradRev, as_value_column_array_or_scalar_ref_scalar) {
  using stan::math::as_value_column_array_or_scalar_ref;
  stan::math::var a = 2.5;
  double a_val = as_value_column_array_or_scalar_ref(a);
  EXPECT_EQ(2.5, a_val);
}

TEST_F(AgradRev, as_value_column_array_or_scalar_ref_reads_vari) {
  using stan::math::var;
  Eigen::Matrix<var, -1, 1> a = Eigen::VectorXd::Ones(3);
  Eigen::Matrix<var, 1, -1> b = Eigen::RowVectorXd::Ones(3);
  std::vector<var> c{var(1.0), var(1.0), var(1.0)};
  expect_view_reads_vari(a);
  expect_view_reads_vari(b);
  expect_view_reads_vari(c);
}

TEST_F(AgradRev, as_value_column_array_or_scalar_ref_varmat) {
  using stan::math::as_value_column_array_or_scalar_ref;
  Eigen::VectorXd x(3);
  x << 1.5, -2, 3;
  stan::math::var_value<Eigen::VectorXd> a = x;
  decltype(auto) a_val = as_value_column_array_or_scalar_ref(a);
  EXPECT_MATRIX_EQ(x.array(), a_val);
  a.vi_->val_.coeffRef(1) = 7;
  EXPECT_FLOAT_EQ(7, a_val[1]);
}
