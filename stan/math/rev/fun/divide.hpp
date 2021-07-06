#ifndef STAN_MATH_REV_FUN_DIVIDE_HPP
#define STAN_MATH_REV_FUN_DIVIDE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/fun/to_arena.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {
/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat A type inheriting from `EigenBase` with an `Arithmetic` scalar
 * type.
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, require_eigen_vt<std::is_arithmetic, Mat>* = nullptr>
inline auto divide(const Mat& m, var c) {
  auto inv_c = (1.0 / c.val());
  arena_t<promote_scalar_t<var, Mat>> res = inv_c * m.array();
  reverse_pass_callback([c, inv_c, res]() mutable {
    c.adj() -= inv_c * (res.adj().array() * res.val().array()).sum();
  });
  return promote_scalar_t<var, Mat>(res);
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat Either a type inheriting from `EigenBase` with a scalar type of
 * `var` or a `var_value<T>` with type `T` inheriting from `EigenBase`.
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, require_matrix_st<is_var, Mat>* = nullptr>
inline auto divide(const Mat& m, double c) {
  arena_t<promote_scalar_t<var, Mat>> arena_m = m;
  auto inv_c = (1.0 / c);
  arena_t<promote_scalar_t<var, Mat>> res = inv_c * arena_m.val();
  reverse_pass_callback([c, inv_c, arena_m, res]() mutable {
    arena_m.adj().array() += inv_c * res.adj_op().array();
  });
  return promote_scalar_t<var, Mat>(res);
}

/**
 * Return matrix divided by scalar.
 *
 * @tparam Mat Either a type inheriting from `EigenBase` with a scalar type of
 * `var` or a `var_value<T>` with type `T` inheriting from `EigenBase`.
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, require_matrix_st<is_var, Mat>* = nullptr>
inline auto divide(const Mat& m, var c) {
  arena_t<plain_type_t<Mat>> arena_m = m;
  auto inv_c = (1.0 / c.val());
  arena_t<plain_type_t<Mat>> res = inv_c * arena_m.val();
  reverse_pass_callback([c, inv_c, arena_m, res]() mutable {
    c.adj() -= inv_c * (res.adj().array() * res.val().array()).sum();
    arena_m.adj().array() += inv_c * res.adj().array();
  });
  return plain_type_t<Mat>(res);
}

/**
 * Return scalar divided by matrix.
 *
 * @tparam Mat Either a type inheriting from `EigenBase` with a scalar type of
 * `var` or a `var_value<T>` with type `T` inheriting from `EigenBase`.
 * @param[in] m specified matrix or expression
 * @param[in] c specified scalar
 * @return matrix divided by the scalar
 */
template <typename Mat, require_matrix_st<is_var, Mat>* = nullptr>
inline auto divide(var c, const Mat& m) {
  arena_t<plain_type_t<Mat>> arena_m = m;
  auto inv_m = to_arena(1.0 / arena_m.val().array());
  arena_t<plain_type_t<Mat>> res = c.val() * inv_m;
  reverse_pass_callback([c, inv_m, arena_m, res]() mutable {
    arena_m.adj().array() -= inv_m * res.adj().array() * res.val().array();
    c.adj() += (inv_m * res.adj().array()).sum();
  });
  return plain_type_t<Mat>(res);
}

/**
 * Return scalar divided by matrix.
 *
 * @tparam Mat A type inheriting from `EigenBase` with an `Arithmetic` scalar
 * type.
 * @param[in] c specified scalar
 * @param[in] m specified matrix or expression
 * @return matrix divided by the scalar
 */
template <typename Mat, require_eigen_vt<std::is_arithmetic, Mat>* = nullptr>
inline auto divide(var c, const Mat& m) {
  arena_t<Mat> arena_m = m;
  auto inv_m = to_arena(1.0 / arena_m.val().array());
  arena_t<promote_scalar_t<var, Mat>> res = c.val() * inv_m;
  reverse_pass_callback([c, inv_m, res]() mutable {
    c.adj() += (inv_m * res.adj().array()).sum();
  });
  return promote_scalar_t<var, Mat>(res);
}

/**
 * Return scalar divided by matrix.
 *
 * @tparam Mat Either a type inheriting from `EigenBase` with a scalar type of
 * `var` or a `var_value<T>` with type `T` inheriting from `EigenBase`.
 * @param[in] c specified scalar
 * @param[in] m specified matrix or expression
 * @return matrix divided by the scalar
 */
template <typename Mat, require_matrix_st<is_var, Mat>* = nullptr>
inline auto divide(double c, const Mat& m) {
  arena_t<promote_scalar_t<var, Mat>> arena_m = m;
  auto inv_m = to_arena(1.0 / arena_m.val().array());
  arena_t<promote_scalar_t<var, Mat>> res = c * inv_m;
  reverse_pass_callback([inv_m, arena_m, res]() mutable {
    arena_m.adj().array() -= inv_m * res.adj().array() * res.val().array();
  });
  return promote_scalar_t<var, Mat>(res);
}

///

/**
 * Return a matrix divided by a matrix elementwise.
 * @tparam Mat1 Either a type inheriting from `EigenBase` or a `var_value<T>`
 *  with type `T` inheriting from `EigenBase`.
 * @tparam Mat2 Either a type inheriting from `EigenBase` or a `var_value<T>`
 *  with type `T` inheriting from `EigenBase`.
 * @param[in] m1 specified matrix or expression
 * @param[in] m2 specified matrix or expression
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_st<is_var_or_arithmetic, Mat1, Mat2>* = nullptr,
          require_any_matrix_st<is_var, Mat1, Mat2>* = nullptr>
inline auto divide(const Mat1& m1, const Mat2& m2) {
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    auto inv_m2 = to_arena(arena_m2.val().array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.val().array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.val().array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      auto inv_times_res = (inv_m2 * res.adj().array()).eval();
      arena_m1.adj().array() += inv_times_res;
      arena_m2.adj().array() -= inv_times_res * res.val().array();
    });
    return ret_type(res);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    auto inv_m2 = to_arena(arena_m2.val().array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      arena_m2.adj().array() -= inv_m2 * res.adj().array() * res.val().array();
    });
    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, Mat2>> arena_m2 = value_of(m2);
    auto inv_m2 = to_arena(arena_m2.array().inverse());
    using val_ret = decltype((inv_m2 * arena_m1.val().array()).matrix().eval());
    using ret_type = return_var_matrix_t<val_ret, Mat1, Mat2>;
    arena_t<ret_type> res = (inv_m2.array() * arena_m1.val().array()).matrix();
    reverse_pass_callback([inv_m2, arena_m1, arena_m2, res]() mutable {
      arena_m1.adj().array() += inv_m2 * res.adj().array();
    });
    return ret_type(res);
  }
}
}  // namespace math
}  // namespace stan

#endif
