#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
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

  // TODO(Sean): This behavior doesn't seem to make sense prima facie
  // Likely change this so that we would expect `double` here instead.
  expect_same_type<MatrixXd, scalar_type<vector<MatrixXd> >::type>();
}
