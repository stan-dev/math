#ifndef TEST_UNIT_MATH_AD_TAPE_UTIL_HPP
#define TEST_UNIT_MATH_AD_TAPE_UTIL_HPP

#include <stan/math/rev/core.hpp>
#include <cstddef>

namespace stan {
namespace test {

/**
 * Number of varis on the reverse-mode stack whose `chain()` runs in the
 * reverse pass.
 */
inline std::size_t chain_stack_size() {
  return stan::math::ChainableStack::instance_->var_stack_.size();
}

/**
 * Number of varis on the reverse-mode stack that are skipped by the reverse
 * pass (values without a chain method, e.g. results of a callback).
 */
inline std::size_t nochain_stack_size() {
  return stan::math::ChainableStack::instance_->var_nochain_stack_.size();
}

/**
 * Total number of varis on the reverse-mode stack.
 */
inline std::size_t tape_size() {
  return chain_stack_size() + nochain_stack_size();
}

/**
 * Records the stack sizes at construction and reports how many chain and
 * nochain varis were added since.
 */
struct tape_delta {
  std::size_t chain_start_{chain_stack_size()};
  std::size_t nochain_start_{nochain_stack_size()};
  std::size_t chain() const { return chain_stack_size() - chain_start_; }
  std::size_t nochain() const { return nochain_stack_size() - nochain_start_; }
  std::size_t total() const { return chain() + nochain(); }
};

}  // namespace test
}  // namespace stan

#endif
