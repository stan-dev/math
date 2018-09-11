#ifndef STAN_MATH_PRIM_ARR_FUN_ASC_LIN_SEQ_BY_ARRAY_HPP
#define STAN_MATH_PRIM_ARR_FUN_ASC_LIN_SEQ_BY_ARRAY_HPP

#include <stan/math/prim/scal/err/check_greater.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <vector>
#include <math.h>
#include <iostream>

namespace stan {
namespace math {

/**
 * Return the specified standard vector in ascending order.
 *
 * @tparam T Type of elements contained in vector.
 * @param xs Vector to order.
 * @return Vector in ascending order.
 * @throw std::domain_error If any of the values are NaN.
 */
template <typename result_t, typename T1, typename T2, typename T3>
inline std::vector<result_t> 
asc_lin_seq_by_array_helper(T1 min, T2 max, T3 by, bool inc_max) {
  check_not_nan("asc_lin_seq_by_array", "max", max);
  check_not_nan("asc_lin_seq_by_array", "min", min);
  check_not_nan("asc_lin_seq_by_array", "by", by);
  check_finite("asc_lin_seq_by_array", "min", min); 
  check_finite("asc_lin_seq_by_array", "max", max); 
  check_positive_finite("asc_lin_seq_by_array", "by", by); 
  check_greater("asc_lin_seq_by_array", "max", max, min); 

  std::vector<result_t> v_seq(ceil((max - min) / static_cast<double>(by)));
  for (size_t i = 0; i < length(v_seq); i++) {
    v_seq[i] = min + i * by;
  }
  if (inc_max) {
    v_seq.push_back(max);
  }
  return v_seq;
}

inline std::vector<int> 
asc_lin_seq_by_array(int min, int max, int by, bool inc_max) {
  return asc_lin_seq_by_array_helper<int>(min, max, by, inc_max);
}

inline std::vector<int> 
asc_lin_seq_by_array(int min, int max, int by) {
  return asc_lin_seq_by_array_helper<int>(min, max, by, false);
}

template <typename T1, typename T2, typename T3>
inline std::vector<typename return_type<T1, T2, T3>::type> 
asc_lin_seq_by_array(T1 min, T2 max, T3 by, bool inc_max) {
  typedef typename return_type<T1, T2, T3>::type result_t;
  return asc_lin_seq_by_array_helper<result_t>(min, max, by, inc_max);
}

template <typename T1, typename T2, typename T3>
inline std::vector<typename return_type<T1, T2, T3>::type> 
asc_lin_seq_by_array(T1 min, T2 max, T3 by) {
  typedef typename return_type<T1, T2, T3>::type result_t;
  return asc_lin_seq_by_array_helper<result_t>(min, max, by, false);
}

}  // namespace math
}  // namespace stan
#endif
