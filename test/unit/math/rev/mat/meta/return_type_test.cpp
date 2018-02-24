#include <stan/math/rev/mat.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <vector>

using stan::math::matrix_d;
using stan::math::matrix_v;
using stan::math::var;
using stan::math::vector_d;
using stan::math::vector_v;
using stan::return_type;
using std::vector;

TEST(MetaTraits, ReturnTypeVarMat) {
  test::expect_same_type<var, return_type<vector_v>::type>();
  test::expect_same_type<var, return_type<matrix_v>::type>();
  test::expect_same_type<var, return_type<matrix_v, double>::type>();
  test::expect_same_type<var, return_type<matrix_v, var>::type>();
  test::expect_same_type<var, return_type<matrix_d, matrix_v>::type>();
}

TEST(MetaTraits, ReturnTypeMatMultivar) {
  // test::expect_same_type<var, return_type<vector<vector_v> >::type>();
  test::expect_same_type<var, return_type<vector<matrix_v> >::type>();
  test::expect_same_type<var, return_type<vector<matrix_v>, double>::type>();
  test::expect_same_type<var, return_type<vector<matrix_v>, var>::type>();
  test::expect_same_type<var, return_type<vector<matrix_d>, matrix_v>::type>();
}

TEST(MetaTraits, ReturnTypeDoubleMat) {
  test::expect_same_type<double, return_type<vector_d>::type>();
  test::expect_same_type<double, return_type<matrix_d, double>::type>();
  test::expect_same_type<var, return_type<matrix_d, var>::type>();
}
