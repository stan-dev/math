#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(log_sum_exp_tests, large_values) {
  using stan::math::var;

  Eigen::VectorXd x = Eigen::VectorXd::Random(5);
  stan::math::promote_scalar_t<var, Eigen::VectorXd> xv = x;

  auto y = stan::math::log_softmax(xv);
  auto y_ref = stan::math::log(stan::math::softmax(xv));

  for(size_t i = 0; i < 1/*y.size()*/; ++i) {
    stan::math::set_zero_all_adjoints();
    y(i).grad();
    auto xadj1 = xv.adj().eval();

    stan::math::set_zero_all_adjoints();
    y_ref(i).grad();
    auto xadj2 = xv.adj().eval();

    EXPECT_FLOAT_EQ(y(i).val(), y_ref(i).val()) <<
      "values y_" << i << " not equal" << std::endl <<
      "values: " << y(i).val() << ", reference: " << y_ref(i).val() << std::endl;

    for(size_t j = 0; j < xadj1.size(); ++j) {
      EXPECT_FLOAT_EQ(xadj1(j), xadj2(j)) <<
	"gradients dy_" << i << "/dx_" << j << " not equal" << std::endl <<
	"gradient: " << xadj1(j) << ", reference: " << xadj2(j) << std::endl;
    }
  }
}
