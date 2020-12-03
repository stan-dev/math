#ifndef STAN_MATH_REV_CORE_CHAINABLE_OBJECT_HPP
#define STAN_MATH_REV_CORE_CHAINABLE_OBJECT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/chainable_alloc.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T>
class chainable_object : public chainable_alloc {
 private:
  T obj_;

 public:
  template <typename S>
  chainable_object(S&& obj) : obj_(std::forward<S>(obj)) {}

  T& get() { return obj_; }
};

template <typename T>
auto make_chainable_ptr(T&& obj) {
  auto ptr = new chainable_object<plain_type_t<T>>(std::forward<T>(obj));
  return ptr;
}

}  // namespace math
}  // namespace stan
#endif
