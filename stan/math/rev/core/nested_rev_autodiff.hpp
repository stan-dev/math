#ifndef STAN_MATH_REV_CORE_nested_rev_autodiff_HPP
#define STAN_MATH_REV_CORE_nested_rev_autodiff_HPP

#include <stan/math/rev/core/recover_memory_nested.hpp>
#include <stan/math/rev/core/set_zero_all_adjoints_nested.hpp>
#include <stan/math/rev/core/start_nested.hpp>

namespace stan {
namespace math {

/**
 * A class following the RAII idiom to start and recover nested autodiff scopes.
 * This is the preferred way to use nested autodiff. Example:
 *
 * var a; // allocated normally
 * {
 *    nested_rev_autodiff nested; // Starts nested autodiff
 *
 *    var nested_var; //allocated on the nested stack
 *    // Do stuff on the nested stack
 *
 *    // Nested stack is automatically recovered at the end of scope where
 *    // nested was declared, including exceptions, returns, etc.
 * }
 * var b;
 */
class nested_rev_autodiff {
 public:
  nested_rev_autodiff() { start_nested(); }

  ~nested_rev_autodiff() { recover_memory_nested(); }

  // Prevent undesirable operations
  nested_rev_autodiff(const nested_rev_autodiff&) = delete;
  nested_rev_autodiff& operator=(const nested_rev_autodiff&) = delete;
  void* operator new(std::size_t) = delete;

  /**
   * Reset all adjoint values in this nested stack
   * to zero.
   **/
  inline void set_zero_all_adjoints() { set_zero_all_adjoints_nested(); }
};

}  // namespace math
}  // namespace stan
#endif
