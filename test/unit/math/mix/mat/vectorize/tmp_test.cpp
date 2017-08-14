#include <stan/math/mix/mat.hpp>
#include <stan/math/fwd/mat/vectorize/apply_scalar_binary.hpp>
#include <stan/math/rev/mat/vectorize/apply_scalar_binary.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/vectorize/binary_foo_fun.hpp>
#include <test/unit/math/mix/mat/vectorize/expect_binary_val_deriv_eq.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp> 
#include <vector>

TEST(BinaryTestCase, MixStdVector) {
  using stan::math::fvar;
  using stan::math::var;
  using stan::math::binary_foo;
  using std::vector;
  using stan::math::pow;

  stan::math::recover_memory();
  //fvar<var> x1(0, 0);
  fvar<var> x1(var(0), var(0));
  fvar<var> x2(var(0), var(0));
  vector<int> a1(5, -4);
  vector<int> a2(5, -4);

  //This is applying pow to fvar<var>, int
  fvar<var> result1 = binary_foo(x1, a1[0]); 
  AVEC exp_b = createAVEC(x1.d_);
  VEC exp_bg;
  result1.d_.grad(exp_b, exp_bg);
  std::cout << "Base case 2nd derivative: " << exp_bg[0] << std::endl;
  stan::math::set_zero_all_adjoints();

  /*
    This is applying pow to fvar<var>, vector<int>.
    Commented out line is how it will be called.
    Code being run is essentially the back-end code.
  */
  //vector<fvar<var> > result2 = binary_foo(x2, a2);
  vector<fvar<var> > result2(a2.size());
  //stan::scalar_seq_view<const fvar<var> > x2_vec(x2);
  for (size_t i = 0; i < result2.size(); ++i) {
    result2[i] = binary_foo(x2, a2[i]); 
    //stan::math::recover_memory();
  }
  AVEC exp_t = createAVEC(x2.d_);
  VEC exp_tg;
  std::cout << result2[0].d_ << std::endl;
  std::cout << result2[0].val() << std::endl;
  //result2[0].d_.grad(exp_t, exp_tg);
  result2[0].d_.grad();
  std::cout << "Vector case 2nd derivative: " 
            << x2.d_.adj() << std::endl;
 
  //Commented out line below is how it is tested in the testing framework
  //expect_binary_val_deriv_eq(result1, x1, a1[0], result2[0], x2, a2[0]);
  stan::math::set_zero_all_adjoints();

/*
  stan::math::matrix_v m(2, 2);
  m << 1, 2, 3, 4;
  var a(2);
  stan::math::matrix_v result3 = stan::math::add(m, 2);
  result3(1).grad();
  std::cout << a.adj() << std::endl;
  stan::math::set_zero_all_adjoints();

  var b(2);
  var c(3);
  var d = b + c;
  d.grad();
  std::cout << b.adj() << std::endl;
*/
  std::vector<int> v;
  std::cout << !stan::length(v) << std::endl;
}

