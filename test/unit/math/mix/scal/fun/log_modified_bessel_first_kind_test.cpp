#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>

struct log_modified_bessel_first_kind_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply(const T1& x1, const T2& x2) {
    using stan::math::log_modified_bessel_first_kind;
    return log_modified_bessel_first_kind(x1, x2);
  }
};

TEST(mathMixCore, logModifiedBesselFirstKind) {
  stan::math::test::test_common_args<log_modified_bessel_first_kind_f, false>();
}
