#ifndef STAN_MATH_PRIM_FUN_REP_ARRAY_HPP
#define STAN_MATH_PRIM_FUN_REP_ARRAY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T_ret, typename In, require_std_vector_t<T_ret>* = nullptr>
inline std::vector<plain_type_t<In>> rep_array(const In& x, int n) {
  using T = plain_type_t<In>;
  check_nonnegative("rep_array", "n", n);
  return std::vector<T>(n, x);
}
template <typename In>
inline std::vector<plain_type_t<In>> rep_array(const In& x, int n) {
  return rep_array<std::vector<plain_type_t<In>>>(x, n);
}

template <typename In>
inline std::vector<std::vector<plain_type_t<In>>> rep_array(const In& x, int m,
                                                            int n) {
  using std::vector;
  using T = plain_type_t<In>;
  check_nonnegative("rep_array", "rows", m);
  check_nonnegative("rep_array", "cols", n);
  return vector<vector<T>>(m, vector<T>(n, x));
}

template <typename In>
inline std::vector<std::vector<std::vector<plain_type_t<In>>>> rep_array(
    const In& x, int k, int m, int n) {
  using std::vector;
  using T = plain_type_t<In>;
  check_nonnegative("rep_array", "shelfs", k);
  check_nonnegative("rep_array", "rows", m);
  check_nonnegative("rep_array", "cols", n);
  return vector<vector<vector<T>>>(k, vector<vector<T>>(m, vector<T>(n, x)));
}

}  // namespace math
}  // namespace stan

#endif
