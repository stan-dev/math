#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>
#include <vector>

struct log_modified_bessel_first_kind_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type apply(
      const T1& x1, const T2& x2) {
    using stan::math::log_modified_bessel_first_kind;
    return log_modified_bessel_first_kind(x1, x2);
  }
};

TEST(mathMixCore, logModifiedBesselFirstKind) {
  std::vector<double> xs
      = {12.5, 2.3, 1.2, 0.5, stan::math::positive_infinity()};
  stan::math::test::test_args<log_modified_bessel_first_kind_f, false>(xs);

  // TODO(carpenter):  need to fix boundary conditions in function
  //   50 : tolerance too tight for autodiff, but basically OK
  //   0 : known that won't work for finite diffs test
  //   -0.5  : autodiff doesn't throw
  //   stan::math::not_a_number() : autodiff doesn't throw
  //   stan::math::negative_infinity() : autodiff doesn't throw
}
