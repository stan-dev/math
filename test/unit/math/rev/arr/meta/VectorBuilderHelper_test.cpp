#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(MetaTraits, VectorBuilderHelper_false_true) {
  using std::vector;
  using stan::VectorBuilderHelper;
  using stan::math::var;
  using stan::length;

  std::vector<var> a_std_vector(3);

  VectorBuilderHelper<double,false,true> dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
}

TEST(MetaTraits, VectorBuilderHelper_true_true) {
  using std::vector;
  using stan::VectorBuilderHelper;
  using stan::math::var;
  using stan::length;

  std::vector<var> a_std_vector(3);

  VectorBuilderHelper<double,true,true> dvv2(length(a_std_vector));
  dvv2[0] = 0.0;
  dvv2[1] = 1.0;
  dvv2[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(1.0, dvv2[1]);
  EXPECT_FLOAT_EQ(2.0, dvv2[2]);  
}
