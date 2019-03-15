#include <stan/math/prim/mat.hpp>
#include <stan/math/rev/core/var.hpp>
#include <gtest/gtest.h>
#include <boost/typeof/typeof.hpp>
#include <type_traits>

using Eigen::Matrix;
using stan::math::promote_elements;
using stan::math::var;

TEST(MathFunctionsMatPromote_Elements, doubleMat2doubleMat) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  promote_elements<Matrix<double, 2, 3>, Matrix<double, 2, 3> > p;
  typedef BOOST_TYPEOF(p.promote(m1)) result_t;
  bool same = std::is_same<Matrix<double, 2, 3>, result_t>::value;
  EXPECT_TRUE(same);
}

TEST(MathFunctionsMatPromote_Elements, doubleMat2varMat) {
  stan::math::matrix_d m1(2, 3);
  m1 << 1, 2, 3, 4, 5, 6;
  promote_elements<Matrix<var, 2, 3>, Matrix<double, 2, 3> > p;
  typedef BOOST_TYPEOF(p.promote(m1)) result_t;
  bool same = std::is_same<Matrix<var, 2, 3>, result_t>::value;
  EXPECT_TRUE(same);
}
