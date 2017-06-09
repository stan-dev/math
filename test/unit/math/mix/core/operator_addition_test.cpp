#include <stan/math.hpp>
#include <test/unit/math/mix/mat/util/autodiff_tester.hpp>

struct op_add_f {
  template <typename T1, typename T2>
  static typename boost::math::tools::promote_args<T1, T2>::type
  apply(const T1& x1, const T2& x2) {
    return x1 + x2;
  }
};

TEST(foo, framework) {
  using stan::math::test::test_ad;

  test_ad<op_add_f>(0.5, 1.3, 1.8);
  test_ad<op_add_f>(1, 2.0, 3);
  test_ad<op_add_f>(-1.5, 1.5, 0);


  double pos_inf = stan::math::positive_infinity();
  test_ad<op_add_f>(pos_inf, 0, pos_inf);
  test_ad<op_add_f>(0, pos_inf, pos_inf);
  test_ad<op_add_f>(pos_inf, pos_inf, pos_inf);

  double neg_inf = stan::math::negative_infinity();
  test_ad<op_add_f>(neg_inf, 0, neg_inf);
  test_ad<op_add_f>(0, neg_inf, neg_inf);
  test_ad<op_add_f>(neg_inf, neg_inf, neg_inf);

  double nan = stan::math::not_a_number();
  test_ad<op_add_f>(nan, 1.5, nan);
  test_ad<op_add_f>(0.2, nan, nan);
  test_ad<op_add_f>(nan, nan, nan);
  test_ad<op_add_f>(pos_inf, neg_inf, nan);
  test_ad<op_add_f>(neg_inf, pos_inf, nan);
}


