#ifndef STAN_MATH_PRIM_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_FUN_VALUE_OF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cstddef>
#include <vector>

namespace stan {
namespace math {

/**
 * Inputs that are arithmetic types or containers of airthmetric types
 * are returned from value_of unchanged
 *
 * @tparam T Input type
 * @param[in] x Input argument
 * @return Forwarded input argument
 **/
template <typename T, require_st_arithmetic<T>* = nullptr>
inline T value_of(T&& x) {
  return std::forward<T>(x);
}

template <typename T, require_complex_t<T>* = nullptr,
          require_t<std::is_arithmetic<
              typename std::decay_t<T>::value_type>>* = nullptr>
inline auto value_of(T&& x) {
  return std::forward<T>(x);
}

template <
    typename T, require_complex_t<T>* = nullptr,
    require_not_arithmetic_t<typename std::decay_t<T>::value_type>* = nullptr>
inline auto value_of(T&& x) {
  using inner_t = partials_type_t<typename std::decay_t<T>::value_type>;
  return std::complex<inner_t>{value_of(x.real()), value_of(x.imag())};
}

/**
 * For std::vectors of non-arithmetic types, return a std::vector composed
 * of value_of applied to each element.
 *
 * @tparam T Input element type
 * @param[in] x Input std::vector
 * @return std::vector of values
 **/
template <typename T, require_std_vector_t<T>* = nullptr,
          require_not_st_arithmetic<T>* = nullptr>
inline auto value_of(const T& x) {
  std::vector<plain_type_t<decltype(value_of(std::declval<value_type_t<T>>()))>>
      out;
  out.reserve(x.size());
  for (auto&& x_elem : x) {
    out.emplace_back(value_of(x_elem));
  }
  return out;
}

/**
 * For Eigen matrices and expressions of non-arithmetic types, return an
 *expression that represents the Eigen::Matrix resulting from applying value_of
 *elementwise
 *
 * @tparam EigMat type of the matrix
 *
 * @param[in] M Matrix to be converted
 * @return Matrix of values
 **/
template <typename EigMat, require_eigen_dense_base_t<EigMat>* = nullptr,
          require_not_st_arithmetic<EigMat>* = nullptr>
inline auto value_of(EigMat&& M) {
  return make_holder(
      [](auto& a) {
        return a.unaryExpr([](const auto& scal) { return value_of(scal); });
      },
      std::forward<EigMat>(M));
}

template <typename EigMat, require_eigen_sparse_base_t<EigMat>* = nullptr,
          require_not_st_arithmetic<EigMat>* = nullptr>
inline auto value_of(EigMat&& M) {
  auto&& M_ref = to_ref(M);
  using scalar_t = decltype(value_of(std::declval<value_type_t<EigMat>>()));
  promote_scalar_t<scalar_t, plain_type_t<EigMat>> ret(M_ref.rows(),
                                                       M_ref.cols());
  ret.reserve(M_ref.nonZeros());
  for (int k = 0; k < M_ref.outerSize(); ++k) {
    for (typename std::decay_t<EigMat>::InnerIterator it(M_ref, k); it; ++it) {
      ret.insert(it.row(), it.col()) = value_of(it.valueRef());
    }
  }
  ret.makeCompressed();
  return ret;
}
template <typename EigMat, require_eigen_sparse_base_t<EigMat>* = nullptr,
          require_st_arithmetic<EigMat>* = nullptr>
inline auto value_of(EigMat&& M) {
  return std::forward<EigMat>(M);
}

}  // namespace math
}  // namespace stan

#endif
