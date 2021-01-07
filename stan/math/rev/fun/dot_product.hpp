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
template <typename T1, typename T2, require_all_container_t<T1, T2>* = nullptr,
          require_any_vt_var<T1, T2>* = nullptr>
inline var dot_product(const T1& v1, const T2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<vector_v> v1_arena = as_column_vector_or_scalar(v1);
    arena_t<vector_v> v2_arena = as_column_vector_or_scalar(v2);
    var res(v1_arena.val().dot(v2_arena.val()));
    reverse_pass_callback([v1_arena, v2_arena, res]() mutable {
      for (Eigen::Index i = 0; i < v1_arena.size(); ++i) {
        const auto res_adj = res.adj();
        v1_arena.coeffRef(i).adj() += res_adj * v2_arena.coeffRef(i).val();
        v2_arena.coeffRef(i).adj() += res_adj * v1_arena.coeffRef(i).val();
      }
    });
    return res;
  } else if (!is_constant<T2>::value) {
    arena_t<vector_v> v2_arena = as_column_vector_or_scalar(v2);
    arena_t<Eigen::VectorXd> v1_val_arena
        = value_of(as_column_vector_or_scalar(v1));
    var res(v1_val_arena.dot(v2_arena.val()));
    reverse_pass_callback([v1_val_arena, v2_arena, res]() mutable {
      v2_arena.adj() += res.adj() * v1_val_arena;
    });
    return res;
  } else {
    arena_t<vector_v> v1_arena = as_column_vector_or_scalar(v1);
    arena_t<Eigen::VectorXd> v2_val_arena
        = value_of(as_column_vector_or_scalar(v2));
    var res(v1_arena.val().dot(v2_val_arena));
    reverse_pass_callback([v1_arena, v2_val_arena, res]() mutable {
      v1_arena.adj() += res.adj() * v2_val_arena;
    });
    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
