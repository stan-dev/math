#ifndef STAN_MATH_REV_FUN_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
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
template <typename T1, typename T2, require_all_vector_t<T1, T2>* = nullptr,
          require_not_complex_t<return_type_t<T1, T2>>* = nullptr,
          require_all_not_std_vector_t<T1, T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline var dot_product(T1&& v1, T2&& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);

  if (v1.size() == 0) {
    return 0.0;
  }
  arena_t<T1> v1_arena = std::forward<T1>(v1);
  arena_t<T2> v2_arena = std::forward<T2>(v2);
  if constexpr (is_autodiffable_v<T1, T2>) {
    return make_callback_var(
        v1_arena.val().dot(v2_arena.val()),
        [v1_arena, v2_arena](const auto& vi) mutable {
          const auto res_adj = vi.adj();
          for (Eigen::Index i = 0; i < v1_arena.size(); ++i) {
            v1_arena.adj().coeffRef(i) += res_adj * v2_arena.val().coeff(i);
            v2_arena.adj().coeffRef(i) += res_adj * v1_arena.val().coeff(i);
          }
        });
  } else if constexpr (is_autodiffable_v<T2>) {
    return make_callback_var(v1_arena.dot(v2_arena.val()),
                             [v1_arena, v2_arena](const auto& vi) mutable {
                               v2_arena.adj().array()
                                   += vi.adj() * v1_arena.array();
                             });
  } else {
    return make_callback_var(v1_arena.val().dot(v2_arena.val()),
                             [v1_arena, v2_arena](const auto& vi) mutable {
                               v1_arena.adj().array()
                                   += vi.adj() * v2_arena.val().array();
                             });
  }
}

}  // namespace math
}  // namespace stan
#endif
