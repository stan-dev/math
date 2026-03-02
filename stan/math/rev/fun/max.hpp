#ifndef STAN_MATH_REV_FUN_MAX_HPP
#define STAN_MATH_REV_FUN_MAX_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/fun/is_nan.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/is_nan.hpp>
#include <stan/math/prim/fun/max.hpp>

namespace stan {
namespace math {

/**
 * Returns the maximum value of the two specified scalar arguments, 
 * where at least one argument is a 'var'.
 *
 * @tparam Scal1 type of first argument (must be 'var' if Scal2 is not)
 * @tparam Scal2 type of second argument (must be 'var' if Scal1 is not)
 * @param[in] x first argument
 * @param[in] y second argument
 * @return maximum value of the two arguments
 */
template <typename Scal1, typename Scal2,
          require_any_var_t<Scal1, Scal2>* = nullptr,
          require_all_not_fvar_t<Scal1, Scal2>* = nullptr>
inline var max(Scal1 x, Scal2 y) {
  double x_val = stan::math::value_of(x);
  double y_val = stan::math::value_of(y);
  if (unlikely(is_nan(x_val))) {
    if (unlikely(is_nan(y_val))) {
      return make_callback_var(NOT_A_NUMBER, [x, y](auto& vi) mutable {  
        if constexpr (is_var<Scal1>::value) x.adj() = NOT_A_NUMBER;
        if constexpr (is_var<Scal2>::value) y.adj() = NOT_A_NUMBER;
      });
    }
    return y;
  }
  if (unlikely(is_nan(y_val))) {
    return x;
  }
  double max_val = std::fmax(x_val, y_val);
  return make_callback_var(max_val, [x, y, x_val, y_val](auto& vi) mutable {
    if (x_val > y_val) {
      if constexpr (is_var<Scal1>::value) x.adj() += vi.adj();
    } else if (y_val > x_val) {
      if constexpr (is_var<Scal2>::value) y.adj() += vi.adj();
    } else {
      if constexpr (is_var<Scal1>::value) x.adj() += vi.adj() * 0.5;
      if constexpr (is_var<Scal2>::value) y.adj() += vi.adj() * 0.5;
    }
  });
}
/**
 * Return the maximum value in a container of 'var'.
 * @tparam T Type of container with 'var' scalar type
 * @param[in] x container (Eigen matrix, vector, row vector or std vector)
 * @return maximum value in the container, or -infinity if size is zero
 */
template <typename T,
          require_st_var<T>* = nullptr,
          require_container_t<T>* = nullptr,
          require_not_var_matrix_t<T>* = nullptr>
inline var max(T&& x) {
  if (unlikely(x.size() == 0)) {
    return NEGATIVE_INFTY;
  }
  return apply_vector_unary<T>::apply(std::forward<T>(x), [](auto&& v) {
    const auto& x_ref = to_ref(v);
    arena_t<std::decay_t<decltype(x_ref)>> x_arena(x_ref);
    double max_val = x_arena.val().maxCoeff();
    return make_callback_var(max_val, [x_arena, max_val](auto& vi) mutable {
      if (unlikely(is_nan(max_val))) {
        for (Eigen::Index i = 0; i < x_arena.size(); ++i) {
          if (is_nan(x_arena.val().coeff(i))) {
            x_arena.adj().coeffRef(i) = NOT_A_NUMBER;
          }
        }
        return;
      }
      auto is_max = (x_arena.val().array() == max_val);
      double count = is_max.template cast<double>().sum();
      double adj_to_propagate = vi.adj() / count;
      x_arena.adj().array() += is_max.template cast<double>() * adj_to_propagate;
    });
  });
}

/**
 * Return the maximum value in a `var_value` containing an Eigen matrix,
 * vector, or row vector.
 *
 * @tparam T Type of `var_value` with an Eigen inner type
 * @param[in] x `var_value` matrix
 * @return maximum value in the matrix, or -infinity if size is zero
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
inline var max(const T& x) {
  if (unlikely(x.size() == 0)) {
    return NEGATIVE_INFTY;
  }
  double max_val = x.val().maxCoeff();
  return make_callback_var(max_val, [x, max_val](auto& vi) mutable {
    if (unlikely(is_nan(max_val))) {
      for (Eigen::Index i = 0; i < x.size(); ++i) {
        if (is_nan(x.val().coeff(i))) {
          x.adj().coeffRef(i) = NOT_A_NUMBER;
        }
      }
      return;
    }
    auto is_max = (x.val().array() == max_val);
    double count = is_max.template cast<double>().sum();
    double adj_to_propagate = vi.adj() / count;
    x.adj().array() += is_max.template cast<double>() * adj_to_propagate;
  });
}

}  // namespace math
}  // namespace stan

#endif
