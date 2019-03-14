
#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>










TEST(MetaTraits, ScalarTypeScal) {
  test::expect_same_type<double, stan::scalar_type<double>::type>();
  test::expect_same_type<int, stan::scalar_type<int>::type>();
}





TEST(MetaTraits_arr, ScalarTypeArray) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double, scalar_type<vector<double> >::type>();
  test::expect_same_type<int, scalar_type<vector<int> >::type>();
  test::expect_same_type<double, scalar_type<vector<vector<double> > >::type>();
}

TEST(MetaTraits_arr, ScalarTypeArrayConst) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<double, scalar_type<const vector<double> >::type>();
  test::expect_same_type<int, scalar_type<const vector<int> >::type>();
  test::expect_same_type<double,
                         scalar_type<const vector<vector<double> > >::type>();
}

TEST(MetaTraits_arr, ScalarTypeArrayConstConst) {
  using stan::scalar_type;
  using std::vector;

  test::expect_same_type<const double,
                         scalar_type<const vector<const double> >::type>();
  test::expect_same_type<const int,
                         scalar_type<const vector<const int> >::type>();
  test::expect_same_type<
      const double,
      scalar_type<const vector<const vector<const double> > >::type>();
}





TEST(MetaTraits_mat, ScalarTypeMat) {
  using Eigen::MatrixXd;
  using Eigen::RowVectorXd;
  using Eigen::VectorXd;
  using stan::scalar_type;
  using std::vector;
  using test::expect_same_type;

  expect_same_type<double, scalar_type<MatrixXd>::type>();
  expect_same_type<double, scalar_type<VectorXd>::type>();
  expect_same_type<double, scalar_type<RowVectorXd>::type>();
  expect_same_type<double, scalar_type<vector<double> >::type>();
  expect_same_type<double, scalar_type<vector<MatrixXd> >::type>();
}
