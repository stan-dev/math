#ifndef STAN_MATH_FWD_FUN_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_FWD_FUN_MDIVIDE_LEFT_TRI_LOW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/mdivide_left.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_left_tri_low.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/subtract.hpp>
namespace stan {
namespace math {

template <typename T1, typename T2,
          require_all_eigen_vt<is_fvar, T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, std::decay_t<T1>::RowsAtCompileTime,
                     std::decay_t<T2>::ColsAtCompileTime>
mdivide_left_tri_low(T1&& A, T2&& b) {
  constexpr int S1 = std::decay_t<T1>::RowsAtCompileTime;
  constexpr int C2 = std::decay_t<T2>::ColsAtCompileTime;

  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }
  decltype(auto) b_ref = to_ref(std::forward<T2>(b));
  decltype(auto) A_ref = to_ref(std::forward<T1>(A));
  auto inv_A_mult_b
      = eval(mdivide_left_tri<Eigen::Lower>(A_ref.val(), b_ref.val()));
  return to_fvar(
      inv_A_mult_b,
      subtract(mdivide_left_tri<Eigen::Lower>(A_ref.val(), b_ref.d()),
               multiply(mdivide_left_tri<Eigen::Lower>(
                            A_ref.val(),
                            A_ref.d().template triangularView<Eigen::Lower>()),
                        inv_A_mult_b)));
}

template <typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_vt_same<double, T1>* = nullptr,
          require_eigen_vt<is_fvar, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T2>, std::decay_t<T1>::RowsAtCompileTime,
                     std::decay_t<T2>::ColsAtCompileTime>
mdivide_left_tri_low(T1&& A, T2&& b) {
  constexpr int S1 = std::decay_t<T1>::RowsAtCompileTime;
  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }
  decltype(auto) A_ref = to_ref(std::forward<T1>(A));
  decltype(auto) b_ref = to_ref(std::forward<T2>(b));
  return to_fvar(mdivide_left_tri<Eigen::Lower>(A_ref, b_ref.val()),
                 mdivide_left_tri<Eigen::Lower>(A_ref, b_ref.d()));
}

template <typename T1, typename T2, require_eigen_vt<is_fvar, T1>* = nullptr,
          require_eigen_vt<std::is_floating_point, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, std::decay_t<T1>::RowsAtCompileTime,
                     std::decay_t<T2>::ColsAtCompileTime>
mdivide_left_tri_low(T1&& A, T2&& b) {
  constexpr int S1 = std::decay_t<T1>::RowsAtCompileTime;
  constexpr int C2 = std::decay_t<T2>::ColsAtCompileTime;
  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }
  decltype(auto) A_ref = to_ref(std::forward<T1>(A));
  auto inv_A_mult_b
      = eval(mdivide_left_tri<Eigen::Lower>(A_ref.val(), std::forward<T2>(b)));
  return to_fvar(
      inv_A_mult_b,
      -multiply(
          mdivide_left_tri<Eigen::Lower>(
              A_ref.val(), A_ref.d().template triangularView<Eigen::Lower>()),
          inv_A_mult_b));
}

}  // namespace math
}  // namespace stan
#endif
