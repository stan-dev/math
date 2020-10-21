#ifndef STAN_MATH_REV_CORE_NEEDS_DESTRUCTOR_HPP
#define STAN_MATH_REV_CORE_NEEDS_DESTRUCTOR_HPP

namespace stan {
namespace math {

/**
 * A base for any type that can be allocated in arena that is not default
 * destructible. This ensures the destructor is actually called when the arena
 * memory is recovered. Different than `chainable_alloc` this does not require
 * objects to be directly allocated by `new`. They can for example also be just
 * members of a larger object.
 */
class needs_destructor {
  needs_destructor** stack_location_;

 public:
  /**
   * Default constructor. Needs to be called by constructors of derived types
   */
  needs_destructor() {
    if (ChainableStack::instance_->memalloc_.in_stack(this)) {
      ChainableStack::instance_->destructor_stack_.push_back(this);
      stack_location_ = &ChainableStack::instance_->destructor_stack_.back();
    } else {
      stack_location_ = nullptr;
    }
  }

  /**
   * Copy constructor.
   * @param other Object to copy
   */
  needs_destructor(const needs_destructor& other) {
    if (ChainableStack::instance_->memalloc_.in_stack(this)) {
      ChainableStack::instance_->destructor_stack_.push_back(this);
      stack_location_ = &ChainableStack::instance_->destructor_stack_.back();
    } else {
      stack_location_ = nullptr;
    }
  }

  /**
   * Destructor. Unschedules the destructor call that would happen when the
   * arena memory is recoverd.
   */
  virtual ~needs_destructor() {
    if (stack_location_ != nullptr) {
      *stack_location_ = nullptr;
    }
  }
};

}  // namespace math
}  // namespace stan
#endif
