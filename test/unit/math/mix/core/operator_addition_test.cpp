#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>

struct op_add_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply(const T1& x1, const T2& x2) {
    return x1 + x2;
  }
};

TEST(mathMixCore, opratorPlus) {
  stan::math::test::test_common_args<op_add_f, false>();
}


