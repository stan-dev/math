#include <stan/math/mix/mat.hpp>
#include <gtest/gtest.h>

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
