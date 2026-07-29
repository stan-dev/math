#ifndef STAN_MATH_REV_FUNCTOR_BROADCAST_ARRAY_HPP
#define STAN_MATH_REV_FUNCTOR_BROADCAST_ARRAY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/sum.hpp>
#include <stan/math/prim/functor/broadcast_array.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename T>
class broadcast_array<T, require_st_var<T>> {
 private:
  std::reference_wrapper<T> prim_;

 public:
  template <typename TT>
  explicit broadcast_array(TT&& prim) : prim_(std::forward<TT>(prim)) {}

  T& operator[](int /*i*/) { return prim_.get(); }

  /** \ingroup type_trait
   * Broadcast array can be assigned a scalar or a vector. If assigned a scalar,
   * it will be used directly. If assigned a vector, the argument will be summed
   * first.
   */
  template <typename Y>
  void operator=(const Y& m) {
    prim_.get() = sum(m);
  }
};

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
