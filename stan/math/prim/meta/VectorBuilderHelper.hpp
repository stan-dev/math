#ifndef STAN_MATH_PRIM_META_VECTORBUILDER_HELPER_HPP
#define STAN_MATH_PRIM_META_VECTORBUILDER_HELPER_HPP

#include <stan/math/prim/meta/contains_vector.hpp>
#include <stdexcept>
#include <stan/math/prim/meta/VectorBuilderHelper.hpp>
#include <vector>





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
class VectorBuilderHelper {
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
class VectorBuilderHelper<T1, true, false> {
 private:
  T1 x_;

 public:
  explicit VectorBuilderHelper(size_t /* n */) : x_(0) {}
  T1& operator[](size_t /* i */) { return x_; }

  typedef T1 type;

  inline type& data() { return x_; }
};

}  // namespace stan








namespace stan {

/**
 * Template specialization for using a vector
 */
template <typename T1>
class VectorBuilderHelper<T1, true, true> {
 private:
  std::vector<T1> x_;

 public:
  explicit VectorBuilderHelper(size_t n) : x_(n) {}

  typedef std::vector<T1> type;

  T1& operator[](size_t i) { return x_[i]; }

  inline type& data() { return x_; }
};
}  // namespace stan
#endif
