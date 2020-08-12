#ifndef STAN_MATH_REV_FUN_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <type_traits>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the dot product.
 *
 * @tparam T1 type of elements in the first vector
 * @tparam T2 type of elements in the second vector
 *
 * @param[in] v1 First vector.
 * @param[in] v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error if sizes of v1 and v2 do not match.
 */
template <typename T1, typename T2, require_all_container_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr>
inline return_type_t<T1, T2> dot_product(const T1& v1, const T2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);

  const auto& v1_col = as_column_vector_or_scalar(v1);
  const auto& v2_col = as_column_vector_or_scalar(v2);

  const auto& v1_val = to_arena_if<!is_constant<T2>::value>(value_of(v1_col));
  const auto& v2_val = to_arena_if<!is_constant<T1>::value>(value_of(v2_col));

  var res(new vari(dot_product(v1_val, v2_val)));

  auto v1_arena = to_arena_if<!is_constant<T1>::value>(v1_col);
  auto v2_arena = to_arena_if<!is_constant<T2>::value>(v2_col);

  reverse_pass_callback([=]() mutable {
    if (!is_constant<T1>::value) {
      forward_as<vector_v>(v1_arena).adj() += res.adj() * v2_val;
    }
    if (!is_constant<T2>::value) {
      forward_as<vector_v>(v2_arena).adj() += res.adj() * v1_val;
    }
  });

  return res;
}

}  // namespace math
}  // namespace stan
#endif
