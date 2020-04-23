#include <test/unit/math/test_ad.hpp>
#include <limits>

template <typename T>
void expect_isinf() {
  // C++ idiom for clients of the math lib
  // std::isnan explicit using; stan::math::nan by ADL
  using std::isinf;
  using std::numeric_limits;
  T inf = numeric_limits<double>::infinity();
  EXPECT_TRUE(isinf(inf));
  EXPECT_TRUE(isinf(-inf));
  T nan = std::numeric_limits<double>::quiet_NaN();
  EXPECT_FALSE(isinf(nan));
  EXPECT_FALSE(isinf(T(1)));
  EXPECT_FALSE(isinf(T(1.0)));
  EXPECT_FALSE(isinf(T(0)));
  EXPECT_FALSE(isinf(T(0.0)));
  EXPECT_FALSE(isinf(T(-1)));
  EXPECT_FALSE(isinf(T(-1.0)));
}

TEST(mixFun, isinf) {
  expect_isinf<d_t>();
  expect_isinf<v_t>();
  expect_isinf<fd_t>();
  expect_isinf<ffd_t>();
  expect_isinf<fv_t>();
  expect_isinf<ffv_t>();
}
