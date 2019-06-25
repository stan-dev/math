
#include <stan/math/prim/meta.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

#include <vector>






using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using stan::return_type;
using std::vector;

TEST(MetaTraits_mat, ReturnTypeMatrixXd) {
  test::expect_same_type<double, return_type<MatrixXd>::type>();
}

TEST(MetaTraits_mat, ReturnTypeVectorXd) {
  test::expect_same_type<double, return_type<VectorXd>::type>();
}

TEST(MetaTraits_mat, ReturnTypeRowVectorXd) {
  test::expect_same_type<double, return_type<RowVectorXd>::type>();
}

TEST(MetaTraits_mat, ReturnTypeArray) {
  test::expect_same_type<double, return_type<vector<int> >::type>();
}
