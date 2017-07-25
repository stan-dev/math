#include <stan/math/mix/mat.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_binary.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_binary.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/vectorize/binary_foo_fun.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints.hpp>
#include <vector>

TEST(BinaryTestCase, MixStdVector) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::binary_foo;
  using std::vector;
  using stan::math::pow;

  fvar<var> x1(0, 0);
  fvar<var> x2(0, 0);
  vector<int> a1(5, -4);
  vector<int> a2(5, -4);

  //This is applying pow to fvar<var>, int
  fvar<var> result1 = binary_foo(x1, a1[0]); 
  AVEC exp_b = createAVEC(x1.d_);
  VEC exp_bg;
  result1.d_.grad(exp_b, exp_bg);
  std::cout << "Base case 2nd derivative: " << exp_bg[0] << std::endl;

  /*
    This is applying pow to fvar<var>, vector<int>.
    Commented out line is how it will be called.
    Code being run is essentially the back-end code.
  */
  //vector<fvar<var> > result2 = binary_foo(x2, a2);
  vector<fvar<var> > result2(a2.size());
  for (size_t i = 0; i < result2.size(); ++i) {
    stan::math::set_zero_all_adjoints();

    /*
      If the line below is commented out, then 2nd derivative of all 
      but the first is nan
    */
    x2 = fvar<var>(0, 0);
    result2[i] = binary_foo(x2, a2[0]); 
    AVEC exp_t = createAVEC(x2.d_);
    VEC exp_tg;
    result2[i].d_.grad(exp_t, exp_tg);
    std::cout << "Vector case " << i << " 2nd derivative: " 
              << exp_tg[0] << std::endl;
  }
 
  //Commented out line below is how it is tested in the testing framework
  //expect_binary_val_deriv_eq(result1, x1, a1[0], result2[0], x2, a2[0]);

}

