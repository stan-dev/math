#ifndef STAN_MATH_PRIM_SCAL_META_VECTORBUILDER_HELPER_HPP
#define STAN_MATH_PRIM_SCAL_META_VECTORBUILDER_HELPER_HPP

#include <stdexcept>

namespace stan {

/**
 *  VectorBuilder allocates type T1 values to be used as
 *  intermediate values. There are 2 template parameters:
 *  - used: boolean variable indicating whether this instance
 *      is used. If this is false, there is no storage allocated
 *      and operator[] throws.
 *  - is_vec: boolean variable indicating whether this instance
 *      should allocate a vector, if it is used. If this is false,
 *      the instance will only allocate a single double value.
 *      If this is true, it will allocate the number requested.
 *      Note that this is calculated based on template parameters
 *      T2 through T7.
 *
 *  These values are mutable.
 */
template <typename T1, bool used, bool is_vec>
class VectorBuilderHelper {  // Base class fails on accessors
 public:
  explicit VectorBuilderHelper(size_t /* n */) {}

  T1& operator[](size_t /* i */) {
    throw std::logic_error("used is false. this should never be called");
  }

  typedef T1 type;

  inline type& data() {
    throw std::logic_error("used is false. this should never be called");
  }
};

template <typename T1>
class VectorBuilderHelper<T1, true, false> {  // When it's used but not a vector
 private:
  T1 x_;

 public:
  explicit VectorBuilderHelper(size_t /* n */) : x_(0) {}
  T1& operator[](size_t /* i */) { return x_; }
  const T1& operator[](size_t /* i */) const { return x_; }
  typedef T1 type;

  inline type& data() { return x_; }
  inline const type& data() const { return x_; }
};

}  // namespace stan
#endif
