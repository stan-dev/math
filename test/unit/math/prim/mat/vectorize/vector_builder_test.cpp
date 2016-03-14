#include <stan/math/prim/mat.hpp>
#include <test/unit/math/prim/mat/vectorize/vector_builder.hpp>
#include <gtest/gtest.h>
#include <vector>


TEST(mathVectorBuilder, test1) {
  using test::math::vector_builder;
  vector_builder<double> x;

  std::vector<double> z = x.build();
  EXPECT_EQ(0, z.size());

  std::vector<double> v = x.add(1).add(2).build();
  EXPECT_EQ(2, v.size());
  EXPECT_FLOAT_EQ(1, v[0]);
  EXPECT_FLOAT_EQ(2, v[1]);
}
