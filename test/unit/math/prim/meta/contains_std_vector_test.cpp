#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMetaPrim, contains_vector_true) {
  using stan::contains_std_vector;
  using std::vector;

  EXPECT_TRUE(contains_std_vector<std::vector<double> >::value);
  EXPECT_TRUE(contains_std_vector<std::vector<int> >::value);
  EXPECT_TRUE(contains_std_vector<const std::vector<double> >::value);
  EXPECT_TRUE(contains_std_vector<const std::vector<int> >::value);
  EXPECT_TRUE(contains_std_vector<std::vector<Eigen::VectorXd> >::value);
  EXPECT_TRUE(contains_std_vector<std::vector<Eigen::RowVectorXd> >::value);
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
  r = contains_std_vector<std::vector<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > >::value;
  EXPECT_TRUE(r);
}

TEST(MathMetaPrim, contains_std_vector_false) {
  using stan::contains_std_vector;

  EXPECT_FALSE(contains_std_vector<double>::value);
  EXPECT_FALSE(contains_std_vector<int>::value);
  EXPECT_FALSE(contains_std_vector<const double>::value);
  EXPECT_FALSE(contains_std_vector<const int>::value);

  EXPECT_FALSE(contains_std_vector<Eigen::VectorXd>::value);
  EXPECT_FALSE(contains_std_vector<Eigen::RowVectorXd>::value);
  bool r = contains_std_vector<
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >::value;
  EXPECT_FALSE(r);
}
