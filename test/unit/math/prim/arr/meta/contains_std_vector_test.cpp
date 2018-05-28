#include <stan/math/prim/arr.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, contains_vector_true) {
  using stan::contains_std_vector;
  using std::vector;

  EXPECT_TRUE(contains_std_vector<std::vector<double> >::value);
  EXPECT_TRUE(contains_std_vector<std::vector<int> >::value);
  EXPECT_TRUE(contains_std_vector<std::vector<const double> >::value);
  EXPECT_TRUE(contains_std_vector<std::vector<const int> >::value);

  bool r = contains_std_vector<std::vector<double>, double, double, double,
                               double, double>::value;
  EXPECT_TRUE(r);
  r = contains_std_vector<double, std::vector<double>, double, double, double,
                          double>::value;
  EXPECT_TRUE(r);
  r = contains_std_vector<double, double, std::vector<double>, double, double,
                          double>::value;
  EXPECT_TRUE(r);
  r = contains_std_vector<double, double, double, std::vector<double>, double,
                          double>::value;
  EXPECT_TRUE(r);
  r = contains_std_vector<double, double, double, double, std::vector<double>,
                          double>::value;
  EXPECT_TRUE(r);
  r = contains_std_vector<double, double, double, double, double,
                          std::vector<double> >::value;
  EXPECT_TRUE(r);
}

TEST(MetaTraits, contains_std_vector_false) {
  using stan::contains_std_vector;

  EXPECT_FALSE(contains_std_vector<double>::value);
  EXPECT_FALSE(contains_std_vector<int>::value);
  EXPECT_FALSE(contains_std_vector<const double>::value);
  EXPECT_FALSE(contains_std_vector<const int>::value);
}
