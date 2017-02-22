#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

using stan::return_type;
using std::vector;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

TEST(MetaTraits, ReturnTypeMatrixXd) {
  test::expect_same_type<double, return_type<MatrixXd>::type>();
}

TEST(MetaTraits, ReturnTypeVectorXd) {
  test::expect_same_type<double, return_type<VectorXd>::type>();
}

TEST(MetaTraits, ReturnTypeRowVectorXd) {
  test::expect_same_type<double, return_type<RowVectorXd>::type>();
}

TEST(MetaTraits, ReturnTypeArray) {
  test::expect_same_type<double, return_type<vector<int> >::type>();
}
