#include <stan/math.hpp>
#include <stan/math/fwd.hpp>
#include <stan/math/mix/fun/typedefs.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>


TEST(AgradFwd, vecExp) {
  using stan::math::exp;
  using stan::math::vector_ffv;
  vector_ffv in1(Eigen::VectorXd::Random(100));
  vector_ffv out_par(100), out_ser(100);
  out_par = exp(in1);
  out_ser = in1.unaryExpr([&](const auto& x){ return exp(x); });
  EXPECT_MATRIX_EQ(out_par.val().val().val(), out_ser.val().val().val());
  EXPECT_MATRIX_EQ(out_par.d().val().val(), out_ser.d().val().val());
}