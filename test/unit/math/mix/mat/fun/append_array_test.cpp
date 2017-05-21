#include <stan/math/mix/mat.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <gtest/gtest.h>

using namespace Eigen;

TEST(AgradMix, append_array_fvar_var_double) {
  std::vector<double> x(2);
  std::vector<stan::math::fvar<stan::math::var > > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(5.0, result[0].val().val());
  EXPECT_FLOAT_EQ(6.0, result[1].val().val());
  EXPECT_FLOAT_EQ(7.0, result[2].val().val());
  EXPECT_FLOAT_EQ(1.0, result[3].val().val());
  EXPECT_FLOAT_EQ(2.0, result[4].val().val());

  EXPECT_FLOAT_EQ(1.5, result[0].tangent().val());
  EXPECT_FLOAT_EQ(-2.5, result[1].tangent().val());
  EXPECT_FLOAT_EQ(-3.5, result[2].tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[3].tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4].tangent().val());

  stan::math::set_zero_all_adjoints();
  result[0].val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().adj());

  stan::math::set_zero_all_adjoints();
  result[0].tangent().grad();
  EXPECT_FLOAT_EQ(0.0, y[1].tangent().adj());
}

TEST(AgradMix, append_array_double_fvar_var) {
  std::vector<double> x(2);
  std::vector<stan::math::fvar<stan::math::var> > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val());

  EXPECT_FLOAT_EQ(0.0, result[0].tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[1].tangent().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val());

  stan::math::set_zero_all_adjoints();
  result[3].val().grad();
  EXPECT_FLOAT_EQ(1.0, y[1].val().adj());

  stan::math::set_zero_all_adjoints();
  result[0].tangent().grad();
  EXPECT_FLOAT_EQ(0.0, y[0].tangent().adj());
}

TEST(AgradMix, append_array_fvar_var_fvar_var) {
  std::vector<stan::math::fvar<stan::math::var> > x(2), y(3), result;

  x[0] = 1.0;
  x[0].d_ = 2.5;
  x[1] = 2.0;
  x[1].d_ = 3.5;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val());

  EXPECT_FLOAT_EQ(2.5, result[0].tangent().val());
  EXPECT_FLOAT_EQ(3.5, result[1].tangent().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val());

  stan::math::set_zero_all_adjoints();
  result[2].val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().adj());

  stan::math::set_zero_all_adjoints();
  result[1].tangent().grad();
  EXPECT_FLOAT_EQ(1.0, x[1].tangent().adj());
}

TEST(AgradMix, append_array_fvar_fvar_var_double) {
  std::vector<double> x(2);
  std::vector<stan::math::fvar<stan::math::fvar<stan::math::var> > > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(5.0, result[0].val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[1].val().val().val());
  EXPECT_FLOAT_EQ(7.0, result[2].val().val().val());
  EXPECT_FLOAT_EQ(1.0, result[3].val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[4].val().val().val());

  EXPECT_FLOAT_EQ(1.5, result[0].tangent().val().val());
  EXPECT_FLOAT_EQ(-2.5, result[1].tangent().val().val());
  EXPECT_FLOAT_EQ(-3.5, result[2].tangent().val().val());
  EXPECT_FLOAT_EQ(0.0, result[3].tangent().val().val());
  EXPECT_FLOAT_EQ(0.0, result[4].tangent().val().val());

  for(size_t i = 0; i < result.size(); i++) {
    EXPECT_FLOAT_EQ(0.0, result[i].val().tangent().val());
    EXPECT_FLOAT_EQ(0.0, result[i].tangent().tangent().val());
  }

  stan::math::set_zero_all_adjoints();
  result[0].val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[0].val().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().tangent().adj());
}

TEST(AgradMix, append_array_double_fvar_fvar_var) {
  std::vector<double> x(2);
  std::vector<stan::math::fvar<stan::math::fvar<stan::math::var> > > y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val().val());

  EXPECT_FLOAT_EQ(0.0, result[0].tangent().val().val());
  EXPECT_FLOAT_EQ(0.0, result[1].tangent().val().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val().val());

  for(size_t i = 0; i < result.size(); i++) {
    EXPECT_FLOAT_EQ(0.0, result[i].val().tangent().val());
    EXPECT_FLOAT_EQ(0.0, result[i].tangent().tangent().val());
  }

  stan::math::set_zero_all_adjoints();
  result[2].val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[2].val().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].val().tangent().adj());
}

TEST(AgradMix, append_array_fvar_fvar_var_fvar_fvar_var) {
  std::vector<stan::math::fvar<stan::math::fvar<stan::math::var> > > x(2), y(3), result;

  x[0] = 1.0;
  x[0].d_ = 2.5;
  x[1] = 2.0;
  x[1].d_ = 3.5;
  y[0] = 5.0;
  y[0].d_ = 1.5;
  y[1] = 6.0;
  y[1].d_ = -2.5;
  y[2] = 7.0;
  y[2].d_ = -3.5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[1].val().val().val());
  EXPECT_FLOAT_EQ(5.0, result[2].val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[3].val().val().val());
  EXPECT_FLOAT_EQ(7.0, result[4].val().val().val());

  EXPECT_FLOAT_EQ(2.5, result[0].tangent().val().val());
  EXPECT_FLOAT_EQ(3.5, result[1].tangent().val().val());
  EXPECT_FLOAT_EQ(1.5, result[2].tangent().val().val());
  EXPECT_FLOAT_EQ(-2.5, result[3].tangent().val().val());
  EXPECT_FLOAT_EQ(-3.5, result[4].tangent().val().val());

  for(size_t i = 0; i < result.size(); i++) {
    EXPECT_FLOAT_EQ(0.0, result[i].tangent().tangent().val());
    EXPECT_FLOAT_EQ(0.0, result[i].val().tangent().val());
  }

  stan::math::set_zero_all_adjoints();
  result[0].tangent().val().grad();
  EXPECT_FLOAT_EQ(1.0, x[0].tangent().val().adj());

  stan::math::set_zero_all_adjoints();
  result[2].tangent().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0].tangent().tangent().adj());
}

TEST(AgradMix, append_array_matrix_double_matrix_fvar_var) {
  std::vector<Matrix<double, Dynamic, Dynamic> > x;
  std::vector<Matrix<stan::math::fvar<stan::math::var>, Dynamic, Dynamic> > y, result;

  for(int i = 0; i < 3; i++)
    x.push_back(Matrix<double, Dynamic, Dynamic>::Zero(3, 3));
  for(int i = 0; i < 2; i++)
    y.push_back(Matrix<stan::math::fvar<stan::math::var>, Dynamic, Dynamic>::Zero(3, 3));

  x[0](0, 0) = 1.0;
  y[1](2, 1) = 2.0;
  y[1](2, 1).d_ = 3.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0, 0).val().val());
  EXPECT_FLOAT_EQ(2.0, result[4](2, 1).val().val());
  EXPECT_FLOAT_EQ(3.0, result[4](2, 1).tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4](2, 2).val().val());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).val().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).tangent().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).tangent().adj());

  for(int i = 0; i < 2; i++)
    y[i] = Matrix<stan::math::fvar<stan::math::var>, Dynamic, Dynamic>::Zero(2, 1);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradMix, append_array_matrix_double_matrix_fvar_fvar_var) {
  std::vector<Matrix<double, Dynamic, Dynamic> > x;
  std::vector<Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic> > y, result;

  for(int i = 0; i < 3; i++)
    x.push_back(Matrix<double, Dynamic, Dynamic>::Zero(3, 3));
  for(int i = 0; i < 2; i++)
    y.push_back(Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic>::Zero(3, 3));

  x[0](0, 0) = 1.0;
  y[1](2, 1).val_ = 2.0;
  y[1](2, 1).val_.d_ = 3.0;
  y[1](2, 1).d_.val_ = 4.0;
  y[1](2, 1).d_.d_ = 5.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0, 0).val().val().val());
  EXPECT_FLOAT_EQ(2.0, result[4](2, 1).val().val().val());
  EXPECT_FLOAT_EQ(3.0, result[4](2, 1).val().tangent().val());
  EXPECT_FLOAT_EQ(4.0, result[4](2, 1).tangent().val().val());
  EXPECT_FLOAT_EQ(5.0, result[4](2, 1).tangent().tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4](2, 2).val().val().val());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).tangent().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).tangent().val().adj());

  for(int i = 0; i < 2; i++)
    y[i] = Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic>::Zero(1, 5);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradMix, append_array_matrix_fvar_fvar_var_matrix_fvar_fvar_var) {
  std::vector<Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic> > x, y, result;

  for(int i = 0; i < 3; i++)
    x.push_back(Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic>::Zero(3, 3));
  for(int i = 0; i < 2; i++)
    y.push_back(Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic>::Zero(3, 3));

  x[0](0, 0).val_ = 1.0;
  x[0](0, 0).val_.d_ = 6.0;
  x[0](0, 0).d_.val_ = 7.0;
  x[0](0, 0).d_.d_ = 8.0;
  y[1](2, 1).val_ = 2.0;
  y[1](2, 1).val_.d_ = 3.0;
  y[1](2, 1).d_.val_ = 4.0;
  y[1](2, 1).d_.d_ = 5.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0, 0).val().val().val());
  EXPECT_FLOAT_EQ(6.0, result[0](0, 0).val().tangent().val());
  EXPECT_FLOAT_EQ(7.0, result[0](0, 0).tangent().val().val());
  EXPECT_FLOAT_EQ(8.0, result[0](0, 0).tangent().tangent().val());
  EXPECT_FLOAT_EQ(2.0, result[4](2, 1).val().val().val());
  EXPECT_FLOAT_EQ(3.0, result[4](2, 1).val().tangent().val());
  EXPECT_FLOAT_EQ(4.0, result[4](2, 1).tangent().val().val());
  EXPECT_FLOAT_EQ(5.0, result[4](2, 1).tangent().tangent().val());
  EXPECT_FLOAT_EQ(0.0, result[4](2, 2).val().val().val());

  stan::math::set_zero_all_adjoints();
  result[0](0, 0).val().val().grad();
  EXPECT_FLOAT_EQ(1.0, x[0](0, 0).val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[0](0, 0).val().tangent().grad();
  EXPECT_FLOAT_EQ(1.0, x[0](0, 0).val().tangent().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).val().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).val().val().adj());

  stan::math::set_zero_all_adjoints();
  result[3](2, 1).tangent().val().grad();
  EXPECT_FLOAT_EQ(1.0, y[0](2, 1).tangent().val().adj());

  for(int i = 0; i < 2; i++)
    y[i] = Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic>::Zero(1, 5);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradMix, append_array_matrix_types) {
  std::vector<Matrix<double, Dynamic, Dynamic> > xddd(3);
  std::vector<Matrix<stan::math::fvar<stan::math::var>, Dynamic, Dynamic> > xfvdd(4), rfvdd;
  std::vector<Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, Dynamic> > xffvdd(5), rffvdd;

  EXPECT_NO_THROW(rfvdd = stan::math::append_array(xddd, xfvdd));
  EXPECT_NO_THROW(rfvdd = stan::math::append_array(xfvdd, xddd));
  EXPECT_NO_THROW(rfvdd = stan::math::append_array(xfvdd, xfvdd));
  EXPECT_NO_THROW(rffvdd = stan::math::append_array(xddd, xffvdd));
  EXPECT_NO_THROW(rffvdd = stan::math::append_array(xffvdd, xddd));
  EXPECT_NO_THROW(rffvdd = stan::math::append_array(xffvdd, xffvdd));

  std::vector<Matrix<double, Dynamic, 1> > xdd1(3);
  std::vector<Matrix<stan::math::fvar<stan::math::var>, Dynamic, 1> > xfvd1(4), rfvd1;
  std::vector<Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, Dynamic, 1> > xffvd1(5), rffvd1;

  EXPECT_NO_THROW(rfvd1 = stan::math::append_array(xdd1, xfvd1));
  EXPECT_NO_THROW(rfvd1 = stan::math::append_array(xfvd1, xdd1));
  EXPECT_NO_THROW(rfvd1 = stan::math::append_array(xfvd1, xfvd1));
  EXPECT_NO_THROW(rffvd1 = stan::math::append_array(xdd1, xffvd1));
  EXPECT_NO_THROW(rffvd1 = stan::math::append_array(xffvd1, xdd1));
  EXPECT_NO_THROW(rffvd1 = stan::math::append_array(xffvd1, xffvd1));

  std::vector<Matrix<double, 1, Dynamic> > xd1d(3);
  std::vector<Matrix<stan::math::fvar<stan::math::var>, 1, Dynamic> > xfv1d(4), rfv1d;
  std::vector<Matrix<stan::math::fvar<stan::math::fvar<stan::math::var> >, 1, Dynamic> > xffv1d(5), rffv1d;

  EXPECT_NO_THROW(rfv1d = stan::math::append_array(xd1d, xfv1d));
  EXPECT_NO_THROW(rfv1d = stan::math::append_array(xfv1d, xd1d));
  EXPECT_NO_THROW(rfv1d = stan::math::append_array(xfv1d, xfv1d));
  EXPECT_NO_THROW(rffv1d = stan::math::append_array(xd1d, xffv1d));
  EXPECT_NO_THROW(rffv1d = stan::math::append_array(xffv1d, xd1d));
  EXPECT_NO_THROW(rffv1d = stan::math::append_array(xffv1d, xffv1d));
}
