#ifndef STAN_MATH_PRIM_META_APPLY_SCALAR_BINARY_HPP
#define STAN_MATH_PRIM_META_APPLY_SCALAR_BINARY_HPP

#include <stan/math/prim/meta/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/meta/is_stan_scalar.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/prim/meta/require_generics.hpp>
#include <vector>

namespace stan {
namespace math {

// Forward declaration to allow specialisations
template <typename T1, typename T2, typename Enable = void,
          typename Enable2 = void>
struct apply_scalar_binary {};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_all_stan_scalar_t<T1, T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return f(x,y);
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_all_eigen_t<T1, T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return x.binaryExpr(y, f).eval();
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_eigen_t<T1>,
                           require_stan_scalar_t<T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return x.unaryExpr([&f,&y](const auto& v){ return f(v, y); }).eval();
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_stan_scalar_t<T1>,
                           require_eigen_t<T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    return y.unaryExpr([&f,&x](const auto& v){ return f(x, v); }).eval();
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2,
                           require_all_std_vector_vt<is_stan_scalar, T1, T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    decltype(auto) x_vec = as_column_vector_or_scalar(x);
    decltype(auto) y_vec = as_column_vector_or_scalar(y);
    using T_return = value_type_t<decltype(x_vec.binaryExpr(y_vec, f))>;
    std::vector<T_return> result(x.size());
    Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
        = x_vec.binaryExpr(y_vec, f);
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_std_vector_vt<is_stan_scalar, T1>,
                           require_stan_scalar_t<T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    decltype(auto) x_vec = as_column_vector_or_scalar(x);
    using T_return = value_type_t<decltype(f(x[0], y))>;
    std::vector<T_return> result(x.size());
    Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
        = x_vec.unaryExpr([&f,&y](const auto& v){ return f(v, y); });
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_stan_scalar_t<T1>,
                           require_std_vector_vt<is_stan_scalar, T2>> {

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    decltype(auto) y_vec = as_column_vector_or_scalar(y);
    using T_return = value_type_t<decltype(f(x, y[0]))>;
    std::vector<T_return> result(y.size());
    Eigen::Map<Eigen::Matrix<T_return, -1, 1>>(result.data(), result.size())
        = y_vec.unaryExpr([&f,&x](const auto& v){ return f(x, v); });
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2,
                           require_all_std_vector_vt<is_container, T1, T2>> {
  using T1_vt = value_type_t<T1>;
  using T2_vt = value_type_t<T2>;

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    using T_return
           = decltype(apply_scalar_binary<T1_vt, T2_vt>::apply(x[0], y[0], f));
    size_t y_size = y.size();
    std::vector<T_return> result(y_size);
    for (size_t i = 0; i < y_size; ++i)
      result[i] = apply_scalar_binary<T1_vt, T2_vt>::apply(x[i], y[i], f);
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_std_vector_vt<is_container, T1>,
                           require_stan_scalar_t<T2>> {
  using T1_vt = value_type_t<T1>;

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    using T_return
            = decltype(apply_scalar_binary<T1_vt, T2>::apply(x[0], y, f));
    size_t x_size = x.size();
    std::vector<T_return> result(x_size);
    for (size_t i = 0; i < x_size; ++i)
      result[i] = apply_scalar_binary<T1_vt, T2>::apply(x[i], y, f);
    return result;
  }

};

template <typename T1, typename T2>
struct apply_scalar_binary<T1, T2, require_stan_scalar_t<T1>,
                           require_std_vector_vt<is_container, T2>> {
  using T2_vt = value_type_t<T2>;

  template <typename F>
  static inline auto apply(const T1& x, const T2& y, const F& f) {
    using T_return
            = decltype(apply_scalar_binary<T1, T2_vt>::apply(x, y[0], f));
    size_t y_size = y.size();
    std::vector<T_return> result(y_size);
    for (size_t i = 0; i < y_size; ++i)
      result[i] = apply_scalar_binary<T1, T2_vt>::apply(x, y[i], f);
    return result;
  }

};

}  // namespace math
}  // namespace stan
#endif
