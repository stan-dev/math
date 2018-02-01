#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, VectorBuilderHelper_false_false) {
  using stan::VectorBuilderHelper;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  VectorBuilderHelper<double, false, false> dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraits, VectorBuilderHelper_true_false) {
  using stan::VectorBuilderHelper;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  VectorBuilderHelper<double, true, false> dvv2(length(a_std_vector));
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(0.0, dvv2[1]);
  EXPECT_FLOAT_EQ(0.0, dvv2[2]);
  double data2(10);
  EXPECT_NO_THROW(data2 = dvv2.data());
  EXPECT_FLOAT_EQ(0.0, data2);
}
