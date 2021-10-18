#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/sin.hpp>
#include <vector>

TEST(AgradRev, multiple_grads) {
  for (int i = 0; i < 100; ++i) {
    stan::math::var a = 2.0;
    stan::math::var b = 3.0 * a;
    stan::math::var c = sin(a) * b;
    // fixes warning regarding unused variable
    c = c;

    stan::math::var nothing;
  }

  stan::math::var d = 2.0;
  stan::math::var e = 3.0;
  stan::math::var f = d * e;

  std::vector<stan::math::var> x{d, e};
  std::vector<double> grad_f;
  f.grad(x, grad_f);

  EXPECT_FLOAT_EQ(3.0, d.adj());
  EXPECT_FLOAT_EQ(2.0, e.adj());

  EXPECT_FLOAT_EQ(3.0, grad_f[0]);
  EXPECT_FLOAT_EQ(2.0, grad_f[1]);
}

TEST(AgradRev, ensure_first_vari_chained) {
  using stan::math::var;

  // Make sure there aren't any varis on stack
  stan::math::recover_memory();

  int N = 10;
  std::vector<var> vars;

  var total = 0.0;
  for (int i = 0; i < N; ++i) {
    vars.push_back(0.0);
    total += vars.back();
  }

  total.grad();

  EXPECT_FLOAT_EQ(0.0, total.val());
  for (int i = 0; i < N; ++i) {
    EXPECT_FLOAT_EQ(1.0, vars[i].adj());
  }
}

namespace stan {
namespace math {

class test_vari : public vari {
 public:
  test_vari() : vari(0.0) {}

  virtual void chain() {
    stan::math::nested_rev_autodiff nested;

    // Add enough vars to make the the var_stack_ vector reallocate
    int N_new_vars = ChainableStack::instance_->var_stack_.capacity() + 1;

    var total = 0.0;
    for (int i = 0; i < N_new_vars; ++i) {
      total += i;
    }

    total.grad();
  }
};

}  // namespace math
}  // namespace stan

TEST(AgradRev, nested_grad_during_chain) {
  using stan::math::var;

  var total = 0.0;
  for (int i = 0; i < 2; ++i) {
    total += i;
  }

  var test_var(new stan::math::test_vari());

  test_var.grad();
}
