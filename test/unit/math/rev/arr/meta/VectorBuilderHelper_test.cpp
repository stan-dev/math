#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, VectorBuilderHelper_false_true) {
  using stan::VectorBuilderHelper;
  using stan::length;
  using stan::math::var;
  using std::vector;

  std::vector<var> a_std_vector(3);

  VectorBuilderHelper<double, false, true> dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraits, VectorBuilderHelper_true_true) {
  using stan::VectorBuilderHelper;
  using stan::length;
  using stan::math::var;
  using std::vector;

  var a_var(1);
  std::vector<var> a_std_vector(3);
  VectorBuilderHelper<double, true, true> dvv1(length(a_var));
  dvv1[0] = 0.0;
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  std::vector<double> data1;
  EXPECT_NO_THROW(data1 = dvv1.data());
  EXPECT_EQ(length(a_var), data1.size());

  VectorBuilderHelper<double, true, true> dvv2(length(a_std_vector));
  dvv2[0] = 0.0;
  dvv2[1] = 1.0;
  dvv2[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv2[0]);
  EXPECT_FLOAT_EQ(1.0, dvv2[1]);
  EXPECT_FLOAT_EQ(2.0, dvv2[2]);
  std::vector<double> data2;
  EXPECT_NO_THROW(data2 = dvv2.data());
  EXPECT_EQ(length(a_std_vector), data2.size());
}
