#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/sin.hpp>

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

TEST(AgradRev, grad_in_reverse_mode) {
  using stan::math::var;

  var total = 0.0;
  for (int i = 0; i < 2; ++i) {
    total += i;
  }

  var test_var(new stan::math::test_vari());

  test_var.grad();
}
