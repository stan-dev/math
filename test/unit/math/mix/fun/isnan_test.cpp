#include <test/unit/math/test_ad.hpp>
#include <limits>

template <typename T>
void expect_isnan() {
  // C++ idiom for clients of the math library
  // std::isnan explicit using; stan::math::nan by ADL
  using std::isnan;
  using std::numeric_limits;
  T inf = numeric_limits<double>::infinity();
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(isnan(inf));
  EXPECT_FALSE(isnan(-inf));
  EXPECT_TRUE(isnan(nan));
  EXPECT_FALSE(isnan(T(1)));
  EXPECT_FALSE(isnan(T(1.0)));
  EXPECT_FALSE(isnan(T(0)));
  EXPECT_FALSE(isnan(T(0.0)));
  EXPECT_FALSE(isnan(T(-1)));
  EXPECT_FALSE(isnan(T(-1.0)));
}

TEST(mixFun, isnan) {
  expect_isnan<d_t>();
  expect_isnan<v_t>();
  expect_isnan<fd_t>();
  expect_isnan<ffd_t>();
  expect_isnan<fv_t>();
  expect_isnan<ffv_t>();
}
