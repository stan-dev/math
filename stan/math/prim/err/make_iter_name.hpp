#ifndef STAN_MATH_PRIM_ERR_MAKE_ITER_NAME_HPP
#define STAN_MATH_PRIM_ERR_MAKE_ITER_NAME_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <string>

namespace stan {
namespace math {
namespace internal {
/**
 * Given a name and index, generate a new name with the index attached.
 * This will produce strings of the form `<name>[<i>]`.
 * @param name The name of the argument.
 * @param i The index.
 */
inline auto make_iter_name(const char* name, size_t i) {
  return std::string(name) + "[" + std::to_string(i + stan::error_index::value)
         + "]";
}
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
