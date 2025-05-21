#ifndef STAN_MATH_REV_FUN_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_DOT_PRODUCT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
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
inline var dot_product(const T1& v1, const T2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);

  if (v1.size() == 0) {
    return 0.0;
  }

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> v1_arena = v1;
    arena_t<promote_scalar_t<var, T2>> v2_arena = v2;
    return make_callback_var(
        v1_arena.val().dot(v2_arena.val()),
        [v1_arena, v2_arena](const auto& vi) mutable {
          const auto res_adj = vi.adj();
          for (Eigen::Index i = 0; i < v1_arena.size(); ++i) {
            v1_arena.adj().coeffRef(i) += res_adj * v2_arena.val().coeff(i);
            v2_arena.adj().coeffRef(i) += res_adj * v1_arena.val().coeff(i);
          }
        });
  } else if (!is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T2>> v2_arena = v2;
    arena_t<promote_scalar_t<double, T1>> v1_val_arena = value_of(v1);
    return make_callback_var(v1_val_arena.dot(v2_arena.val()),
                             [v1_val_arena, v2_arena](const auto& vi) mutable {
                               v2_arena.adj().array()
                                   += vi.adj() * v1_val_arena.array();
                             });
  } else {
    arena_t<promote_scalar_t<var, T1>> v1_arena = v1;
    arena_t<promote_scalar_t<double, T2>> v2_val_arena = value_of(v2);
    return make_callback_var(v1_arena.val().dot(v2_val_arena),
                             [v1_arena, v2_val_arena](const auto& vi) mutable {
                               v1_arena.adj().array()
                                   += vi.adj() * v2_val_arena.array();
                             });
  }
}

}  // namespace math
}  // namespace stan
#endif
