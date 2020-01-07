#include <stan/math/rev/meta.hpp>
#include <stan/math/rev.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <vector>

using stan::math::var;
using stan::return_type;

TEST(MetaTraitsRevScal, ReturnTypeVar) {
  test::expect_same_type<var, return_type<var>::type>();
}

TEST(MetaTraitsRevScal, ReturnTypeVarTenParams) {
  test::expect_same_type<var,
                         return_type<double, var, double, int, double, float,
                                     float, float, var, int>::type>();
}

using std::vector;

TEST(MetaTraitsRevArr, ReturnTypeVarArray) {
  test::expect_same_type<var, return_type<vector<var> >::type>();
  test::expect_same_type<var, return_type<vector<var>, double>::type>();
  test::expect_same_type<var, return_type<vector<var>, double>::type>();
}

TEST(MetaTraitsRevArr, ReturnTypeDoubleArray) {
  test::expect_same_type<double, return_type<vector<double> >::type>();
  test::expect_same_type<double, return_type<vector<double>, double>::type>();
  test::expect_same_type<double, return_type<vector<double>, double>::type>();
}

using stan::math::matrix_d;
using stan::math::matrix_v;
using stan::math::vector_d;
using stan::math::vector_v;

TEST(MetaTraitsRevMat, ReturnTypeVarMat) {
  test::expect_same_type<var, return_type<vector_v>::type>();
  test::expect_same_type<var, return_type<matrix_v>::type>();
  test::expect_same_type<var, return_type<matrix_v, double>::type>();
  test::expect_same_type<var, return_type<matrix_v, var>::type>();
  test::expect_same_type<var, return_type<matrix_d, matrix_v>::type>();
}

TEST(MetaTraitsRevMat, ReturnTypeMatMultivar) {
  // test::expect_same_type<var, return_type<vector<vector_v> >::type>();
  test::expect_same_type<var, return_type<vector<matrix_v> >::type>();
  test::expect_same_type<var, return_type<vector<matrix_v>, double>::type>();
  test::expect_same_type<var, return_type<vector<matrix_v>, var>::type>();
  test::expect_same_type<var, return_type<vector<matrix_d>, matrix_v>::type>();
}

TEST(MetaTraitsRevMat, ReturnTypeDoubleMat) {
  test::expect_same_type<double, return_type<vector_d>::type>();
  test::expect_same_type<double, return_type<matrix_d, double>::type>();
  test::expect_same_type<var, return_type<matrix_d, var>::type>();
}
