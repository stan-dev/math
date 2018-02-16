#include <gtest/gtest.h>

TEST(MathMeta, ensure_c11_features_present) {
  int s[] = {2, 4};
  auto f = [&s](auto i) { return i + s[0]; };
  EXPECT_EQ(4, f(2));
}
