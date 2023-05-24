#ifndef STAN_MATH_PRIM_ERR_MAKE_ITER_NAME_HPP
#define STAN_MATH_PRIM_ERR_MAKE_ITER_NAME_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <string>

namespace stan {
namespace math {
namespace internal {

inline constexpr auto construct_idx() noexcept { return ""; }

template <typename... Idxs>
inline auto construct_idx(size_t i, Idxs... idxs) {
  if (sizeof...(Idxs) > 0) {
    return std::to_string(i + stan::error_index::value) + ", "
           + construct_idx(idxs...);
  } else {
    return std::to_string(i + stan::error_index::value);
  }
}

inline auto make_iter_name(const char* name) { return std::string(name); }

/**
 * Given a name and index, generate a new name with the index attached.
 * This will produce strings of the form `<name>[<i>]`.
 * @param name The name of the argument.
 * @param idxs The index.
 */
template <typename... Idxs>
inline auto make_iter_name(const char* name, Idxs... idxs) {
  return (std::string(name) + "[" + construct_idx(idxs...) + "]");
}
}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
