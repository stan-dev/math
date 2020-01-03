#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <string>
#include <vector>

TEST(std_normal_lccdf, tail) {
  using stan::math::std_normal_lccdf;
  using stan::math::var;
  using std::exp;

  EXPECT_FLOAT_EQ(
      1, -6.661338147750941214694e-16 / (std_normal_lccdf(var(-8.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -3.186340080674249758114e-14 / (std_normal_lccdf(var(-7.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -1.279865102788699562477e-12 / (std_normal_lccdf(var(-7.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -4.015998644826973564545e-11 / (std_normal_lccdf(var(-6.5)).val()));

  EXPECT_FLOAT_EQ(
      1, -9.865877009111571184118e-10 / (std_normal_lccdf(var(-6.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -1.898956265833514866414e-08 / (std_normal_lccdf(var(-5.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -2.866516130081049047962e-07 / (std_normal_lccdf(var(-5.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -3.397678896843115195074e-06 / (std_normal_lccdf(var(-4.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -3.167174337748932124543e-05 / (std_normal_lccdf(var(-4.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -0.0002326561413768195969113 / (std_normal_lccdf(var(-3.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -0.001350809964748202673598 / (std_normal_lccdf(var(-3.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -0.0062290254858600267035 / (std_normal_lccdf(var(-2.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -0.02301290932896348992442 / (std_normal_lccdf(var(-2.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -0.06914345561223400604689 / (std_normal_lccdf(var(-1.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -0.1727537790234499048836 / (std_normal_lccdf(var(-1.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -0.3689464152886565151412 / (std_normal_lccdf(var(-0.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  -0.6931471805599452862268 / (std_normal_lccdf(var(0)).val()));
  EXPECT_FLOAT_EQ(
      1, -1.175911761593618320987 / (std_normal_lccdf(var(0.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -1.841021645009263352222 / (std_normal_lccdf(var(1.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -2.705944400823889317564 / (std_normal_lccdf(var(1.5)).val()));
  EXPECT_FLOAT_EQ(1,
                  -3.78318433368203210776 / (std_normal_lccdf(var(2.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -5.081648277278686620662 / (std_normal_lccdf(var(2.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -6.607726221510342945464 / (std_normal_lccdf(var(3.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -8.366065308344028395027 / (std_normal_lccdf(var(3.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -10.36010148652728979357 / (std_normal_lccdf(var(4.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -12.59241973571053385683 / (std_normal_lccdf(var(4.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -15.06499839383403838156 / (std_normal_lccdf(var(5.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -17.77937635198566113104 / (std_normal_lccdf(var(5.5)).val()));
  EXPECT_FLOAT_EQ(
      1, -20.73676889383495947072 / (std_normal_lccdf(var(6.0)).val()));
  EXPECT_FLOAT_EQ(
      1, -23.93814997800869548428 / (std_normal_lccdf(var(6.5)).val()));

  stan::math::recover_memory();
}

void test_value_and_derivatives(double expected_val, double y_dbl) {
  using stan::math::is_nan;
  using stan::math::std_normal_lccdf;
  using stan::math::var;

  std::stringstream msg_ss;
  msg_ss << "parameters: (" << y_dbl << ")";
  std::string msg = msg_ss.str();
  SCOPED_TRACE(msg);

  var y(y_dbl);
  var val = std_normal_lccdf(y);
  std::vector<var> x;
  x.push_back(y);
  std::vector<double> gradients;
  gradients.clear();
  val.grad(x, gradients);

  double e = 1e-10;
  double inv2e = 0.5 / e;
  double finite_diff
      = (std_normal_lccdf(y_dbl + e) - std_normal_lccdf(y_dbl - e)) * inv2e;

  EXPECT_FLOAT_EQ(expected_val, val.val());
  EXPECT_FALSE(is_nan(gradients[0]));
  if (!is_nan(gradients[0])) {
    if (!is_nan(finite_diff))
      EXPECT_NEAR(finite_diff, gradients[0], 1e-2);
    else
      EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), gradients[0]);
  }
}

TEST(std_normal_lccdf, derivatives) {
  test_value_and_derivatives(log(0.5), 0.0);
  test_value_and_derivatives(0.0, -20.0);
  test_value_and_derivatives(0.0, -30.0);
  test_value_and_derivatives(0.0, -50.0);
  test_value_and_derivatives(-std::numeric_limits<double>::infinity(), 30.0);
}
