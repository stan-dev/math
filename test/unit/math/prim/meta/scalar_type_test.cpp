#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(MathMetaPrim, ScalarTypeScal) {
  test::expect_same_type<double, stan::scalar_type<double>::type>();
  test::expect_same_type<int, stan::scalar_type<int>::type>();
}

TEST(MathMetaPrim, ScalarTypeArray) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double, scalar_type<vector<double>>::type>();
  test::expect_same_type<int, scalar_type<vector<int>>::type>();
  test::expect_same_type<double, scalar_type<vector<vector<double>>>::type>();
}

TEST(MathMetaPrim, ScalarTypeArrayConst) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double, scalar_type<const vector<double>>::type>();
  test::expect_same_type<int, scalar_type<const vector<int>>::type>();
  test::expect_same_type<double,
                         scalar_type<const vector<vector<double>>>::type>();
}

TEST(MathMetaPrim, ScalarTypeArrayConstConst) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double const*,
                         scalar_type<const vector<double const*>>::type>();
  test::expect_same_type<int const*,
                         scalar_type<const vector<int const*>>::type>();
  test::expect_same_type<
      double const*, scalar_type<const vector<vector<double const*>>>::type>();
}

TEST(MathMetaPrim, ScalarTypeMat) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::scalar_type;
  using std::vector;
  using test::expect_same_type;

  expect_same_type<double, scalar_type<MatrixXd>::type>();
  expect_same_type<double, scalar_type<VectorXd>::type>();
  expect_same_type<double, scalar_type<RowVectorXd>::type>();
  expect_same_type<double, scalar_type<vector<double>>::type>();
  expect_same_type<double, scalar_type<vector<MatrixXd>>::type>();
}
