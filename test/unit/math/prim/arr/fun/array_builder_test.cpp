#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::array_builder;
using std::vector;

TEST(MathArray, arrayBuilder) {
  EXPECT_EQ(0U, array_builder<double>().array().size());

  vector<double> v3 = array_builder<double>().add(1).add(3).add(2).array();
  EXPECT_EQ(3U, v3.size());
  EXPECT_FLOAT_EQ(1.0, v3[0]);
  EXPECT_FLOAT_EQ(3.0, v3[1]);
  EXPECT_FLOAT_EQ(2.0, v3[2]);

  vector<vector<int> > v3v2
      = array_builder<vector<int> >()
            .add(array_builder<int>().add(1).add(2).array())
            .add(array_builder<int>().add(3).add(4).array())
            .add(array_builder<int>().add(5).add(6).array())
            .array();

  EXPECT_EQ(3U, v3v2.size());
  for (size_t i = 0; i < 3; ++i)
    EXPECT_EQ(2U, v3v2[i].size());
  EXPECT_EQ(1, v3v2[0][0]);
  EXPECT_EQ(2, v3v2[0][1]);
  EXPECT_EQ(3, v3v2[1][0]);
  EXPECT_EQ(4, v3v2[1][1]);
  EXPECT_EQ(5, v3v2[2][0]);
  EXPECT_EQ(6, v3v2[2][1]);
}
