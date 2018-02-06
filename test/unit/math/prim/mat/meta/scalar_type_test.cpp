#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

TEST(MetaTraits, ScalarTypeMat) {
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
