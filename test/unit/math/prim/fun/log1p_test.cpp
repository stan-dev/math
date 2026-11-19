#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

TEST(MathFunctions, log1p) {
  double x;
  x = 0;
  EXPECT_FLOAT_EQ(0.0, stan::math::log1p(x));
  x = 0.0000001;
  EXPECT_FLOAT_EQ(0.0000001, stan::math::log1p(x));
  x = 0.001;
  EXPECT_FLOAT_EQ(0.0009995003, stan::math::log1p(x));
  x = 0.1;
  EXPECT_FLOAT_EQ(0.09531018, stan::math::log1p(x));
  x = 1;
  EXPECT_FLOAT_EQ(0.6931472, stan::math::log1p(x));
  x = 10;
  EXPECT_FLOAT_EQ(2.397895, stan::math::log1p(x));

  x = -0.0000001;
  EXPECT_FLOAT_EQ(-0.0000001, stan::math::log1p(x));
  x = -0.001;
  EXPECT_FLOAT_EQ(-0.0010005, stan::math::log1p(x));
  x = -0.1;
  EXPECT_FLOAT_EQ(-0.1053605, stan::math::log1p(x));
  x = -0.999;
  EXPECT_FLOAT_EQ(-6.907755, stan::math::log1p(x));
}
TEST(MathFunctions, log1pNegInfReturn) {
  EXPECT_EQ(-std::numeric_limits<double>::infinity(), stan::math::log1p(-1.0));
}

TEST(MathFunctions, log1p_exception) {
  EXPECT_THROW(stan::math::log1p(-10.0), std::domain_error);
}

TEST(MathFunctions, log1p_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_TRUE(std::isnan(stan::math::log1p(nan)));
}

TEST(MathFunctions, log1p_double_containers) {
  using stan::math::log1p;
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> v{nan, -1.0, -0.5, 1e-10, 3.0};
  Eigen::VectorXd ev = Eigen::Map<Eigen::VectorXd>(v.data(), v.size());
  Eigen::RowVectorXd rv = ev.transpose();

  std::vector<double> lv = log1p(v);
  Eigen::VectorXd lev = log1p(ev);
  Eigen::RowVectorXd lrv = log1p(rv);
  Eigen::VectorXd lexpr = log1p(ev * 1.0);
  for (size_t i = 0; i < v.size(); ++i) {
    const double expected = stan::math::log1p(v[i]);
    if (std::isnan(expected)) {
      EXPECT_TRUE(std::isnan(lv[i]));
      EXPECT_TRUE(std::isnan(lev(i)));
      EXPECT_TRUE(std::isnan(lrv(i)));
      EXPECT_TRUE(std::isnan(lexpr(i)));
    } else {
      EXPECT_DOUBLE_EQ(expected, lv[i]);
      EXPECT_DOUBLE_EQ(expected, lev(i));
      EXPECT_DOUBLE_EQ(expected, lrv(i));
      EXPECT_DOUBLE_EQ(expected, lexpr(i));
    }
  }
  EXPECT_EQ(-inf, lev(1));

  std::vector<double> bad{0.5, -2.0};
  EXPECT_THROW(log1p(bad), std::domain_error);
  Eigen::VectorXd bad_e(2);
  bad_e << 0.5, -2.0;
  EXPECT_THROW(Eigen::VectorXd r = log1p(bad_e), std::domain_error);
}
