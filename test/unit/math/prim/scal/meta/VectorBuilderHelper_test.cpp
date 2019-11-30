#include <stan/math/prim/arr.hpp>
#include <stan/math/prim/scal.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraitsPrimScal, VectorBuilderHelper_false_false) {
  using stan::VectorBuilderHelper;
  using stan::length;

  double a_double(1);

  VectorBuilderHelper<double, false, false> dvv1(length(a_double));
  EXPECT_THROW(dvv1[0], std::logic_error);
  EXPECT_THROW(dvv1.data(), std::logic_error);
}

TEST(MetaTraitsPrimScal, VectorBuilderHelper_true_false) {
  using stan::VectorBuilderHelper;
  using stan::length;

  double a_double(1);

  VectorBuilderHelper<double, true, false> dvv1(length(a_double));
  EXPECT_FLOAT_EQ(0.0, dvv1[0]);
  EXPECT_FLOAT_EQ(0.0, dvv1[1]);
  EXPECT_FLOAT_EQ(0.0, dvv1[100]);

  double data = dvv1.data();
  EXPECT_FLOAT_EQ(0.0, data);
}

TEST(MetaTraitsPrimArr, VectorBuilderHelper_false_false) {
  using stan::VectorBuilderHelper;
  using stan::length;
  using std::vector;

  std::vector<double> a_std_vector(3);

  VectorBuilderHelper<double, false, false> dvv2(length(a_std_vector));
  EXPECT_THROW(dvv2[0], std::logic_error);
  EXPECT_THROW(dvv2.data(), std::logic_error);
}

TEST(MetaTraitsPrimArr, VectorBuilderHelper_true_false) {
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
