#include <stan/math/prim.hpp>
#include <test/unit/math/prim/fun/binary_scalar_tester.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathFunctions, log_modified_bessel_first_kind) {
  using stan::math::log_modified_bessel_first_kind;
  using std::log;
  using std::sqrt;

  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(sqrt(3), sqrt(2)),
                  log(boost::math::cyl_bessel_i(sqrt(3), sqrt(2))));

  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(0.5, 10),
                  log(sqrt(2 / (M_PI * 10)) * std::sinh(10)));

  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(-0.5, 10),
                  log(sqrt(2 / (M_PI * 10)) * std::cosh(10)));

  // code branches for z > 100 and v sufficiently small
  EXPECT_FLOAT_EQ(log_modified_bessel_first_kind(.1, 100),
                  log_modified_bessel_first_kind(.1, 100 + 1e-14));

  // code branches when v is -1
  EXPECT_FLOAT_EQ(-0.5706479874908312, log_modified_bessel_first_kind(-1, 1));
  EXPECT_FLOAT_EQ(-0.5706479874908312,
                  log_modified_bessel_first_kind(-1 + 1e-16, 1));

  // code branches at v = 0
  EXPECT_FLOAT_EQ(0.235914358507179,
                  log_modified_bessel_first_kind(0, 1));  // integers promoted
  EXPECT_FLOAT_EQ(0.235914358507179, log_modified_bessel_first_kind(1e-16, 1));

  EXPECT_FLOAT_EQ(5.82456472298118,
                  log_modified_bessel_first_kind(0, 7.75 - 5e-16));
  EXPECT_FLOAT_EQ(5.82456472298118, log_modified_bessel_first_kind(0, 7.75));
  EXPECT_FLOAT_EQ(5.82456472298118,
                  log_modified_bessel_first_kind(1e-16, 7.75));

  EXPECT_FLOAT_EQ(495.974007668107,
                  log_modified_bessel_first_kind(0, 500 - 1e-13));
  EXPECT_FLOAT_EQ(495.974007668107, log_modified_bessel_first_kind(0, 500));
  EXPECT_FLOAT_EQ(495.974007668107, log_modified_bessel_first_kind(1e-16, 500));

  // code branches at v = 1
  EXPECT_FLOAT_EQ(-0.570647987490831,
                  log_modified_bessel_first_kind(1, 1));  // integers promoted
  EXPECT_FLOAT_EQ(-0.570647987490831,
                  log_modified_bessel_first_kind(1 + 1e-16, 1));

  EXPECT_FLOAT_EQ(5.75527527206771,
                  log_modified_bessel_first_kind(1, 7.75 - 5e-16));
  EXPECT_FLOAT_EQ(5.75527527206771, log_modified_bessel_first_kind(1, 7.75));
  EXPECT_FLOAT_EQ(5.75527527206771,
                  log_modified_bessel_first_kind(1 + 1e-16, 7.75));

  EXPECT_FLOAT_EQ(495.973006666268,
                  log_modified_bessel_first_kind(1, 500 - 1e-13));
  EXPECT_FLOAT_EQ(495.973006666268, log_modified_bessel_first_kind(1, 500));
  EXPECT_FLOAT_EQ(495.973006666268,
                  log_modified_bessel_first_kind(1 + 1e-16, 500));

  // limiting cases
  EXPECT_EQ(0, log_modified_bessel_first_kind(0, 0));
  EXPECT_EQ(-std::numeric_limits<double>::infinity(),
            log_modified_bessel_first_kind(1, 0));
  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            log_modified_bessel_first_kind(-0.1, 0));
}

TEST(MathFunctions, log_modified_bessel_first_kind_throw) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  using stan::math::log_modified_bessel_first_kind;

  EXPECT_THROW(log_modified_bessel_first_kind(1, nan), std::domain_error);
  EXPECT_THROW(log_modified_bessel_first_kind(nan, 1), std::domain_error);
  EXPECT_THROW(log_modified_bessel_first_kind(.5, -2), std::domain_error);
  EXPECT_THROW(log_modified_bessel_first_kind(-2, 1), std::domain_error);
}

TEST(MathFunctions, log_modified_bessel_first_kind_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_modified_bessel_first_kind;
    return log_modified_bessel_first_kind(x1, x2);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.8, 3.24, 1.8;
  Eigen::VectorXd in2(3);
  in2 << 7.3, 0.7, 2.8;
  stan::test::binary_scalar_tester(f, in1, in2);
}
