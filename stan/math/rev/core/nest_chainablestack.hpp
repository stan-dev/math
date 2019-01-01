#ifndef STAN_MATH_REV_CORE_NEST_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_NEST_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
#include <stan/math/rev/core/chainablestack.hpp>

#include <stan/math/rev/core/vari.hpp>
#include <stan/math/rev/core/var.hpp>

#include <vector>

namespace stan {
namespace math {

/**
 * A variable implementation which allows to nest a ChainableStack
 * instance as part of another ChainableStack. This is useful if the
 * AD tape is written in independent processes and is then merged
 * in a lazy way by linking the child stacks as part of the parent
 * stack.
 */

class nested_chainablestack_vari : public vari {
 protected:
  ChainableStack::AutodiffStackStorage& nested_;
  const std::size_t start_;
  const std::size_t end_;

 public:
  nested_chainablestack_vari(ChainableStack::AutodiffStackStorage& nested,
                             std::size_t start, std::size_t end)
      : vari(0), nested_(nested), start_(start), end_(end) {}

  void chain() {
    typedef std::vector<vari*>::reverse_iterator it_t;
    // std::cout << "Chaining nested stack " << nested_.id_ << " from " <<
    // start_ << " to " << end_ << std::endl;
    const it_t rend
        = nested_.var_stack_.rbegin() + (nested_.var_stack_.size() - start_);
    for (it_t it
         = nested_.var_stack_.rbegin() + (nested_.var_stack_.size() - end_);
         it != rend; ++it) {
      (*it)->chain();
    }
  }
};

inline void register_nested_chainablestack(
    ChainableStack::AutodiffStackStorage& nested, std::size_t start,
    std::size_t end) {
  new nested_chainablestack_vari(nested, start, end);
}

}  // namespace math
}  // namespace stan

#endif
