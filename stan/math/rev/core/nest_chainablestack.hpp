#ifndef STAN_MATH_REV_CORE_NEST_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_NEST_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/var.hpp>

namespace stan {
namespace math {

/**
 * A variable implementation which allows to nest a ChainableStack
 * instance as part of another ChainableStack. This is useful if the
 * AD tape can be written in independent processes and is then merged
 * in a lazy way by linking the child stacks as part of the parent
 * stack.
 */

class nested_chainablestack_vari : public {
 protected:
  typedef std::vector<vari*>::reverse_iterator it_t;
  const it_t nested_begin_;
  const it_t nested_end_;

 public:
  nested_chainablestack_vari(const it_t nested_begin, const it_t nested_end)
      : nested_begin_(nested_begin), nested_end_(nested_end) {
    // TODO: check that the nested stack is non-nested (or make sure
    // it will work)
  }

  void chain() {
    for (it_t it = nested_begin_; it != nested_end_; ++it) {
      (*it)->chain();
    }
  }
};

inline void register_nested_chainablestack(const it_t nested_begin,
                                           const it_t nested_end) {
  new nested_chainablestack_vari(nested_begin, nested_end);
}

}  // namespace math
}  // namespace stan

#endif
