#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <test/unit/util.hpp>
#include <vector>

TEST(MetaTraits, ScalarTypeMat) {
  using stan::scalar_type;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;
  using Eigen::RowVectorXd;
  using test::expect_same_type;
  using std::vector;

  expect_same_type<double, scalar_type<MatrixXd>::type>();
  expect_same_type<double, scalar_type<VectorXd>::type>();
  expect_same_type<double, scalar_type<RowVectorXd>::type>();
  expect_same_type<double, scalar_type<vector<double> >::type>();
  expect_same_type<double, scalar_type<vector<MatrixXd> >::type>();
}
