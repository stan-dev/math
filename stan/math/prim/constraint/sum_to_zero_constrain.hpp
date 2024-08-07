#ifndef STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP
#define STAN_MATH_PRIM_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined using a modified version of the
 * the inverse of the isometric log ratio transform (ILR).
 * See:
 * Egozcue, Juan Jose; Pawlowsky-Glahn, Vera; Mateu-Figueras, Gloria;
 * Barcelo-Vidal, Carles (2003), "Isometric logratio transformations for
 * compositional data analysis", Mathematical Geology, 35 (3): 279–300,
 * doi:10.1023/A:1023818214614, S2CID 122844634
 *
 * This implementation is closer to the description of the same using "pivot
 * coordinates" in
 * Filzmoser, P., Hron, K., Templ, M. (2018). Geometrical Properties of
 * Compositional Data. In: Applied Compositional Data Analysis. Springer Series
 * in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-319-96422-5_3
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Zero-sum vector of dimensionality K.
 */
template <typename Vec, require_eigen_col_vector_t<Vec>* = nullptr,
          require_not_st_var<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_constrain(const Vec& y) {
  const auto N = y.size();

  plain_type_t<Vec> z = Eigen::VectorXd::Zero(N + 1);
  if (unlikely(N == 0)) {
    return z;
  }

  auto&& y_ref = to_ref(y);

  value_type_t<Vec> sum_w(0);
  for (int i = N; i > 0; --i) {
    double n = static_cast<double>(i);
    auto w = y_ref(i - 1) * inv_sqrt(n * (n + 1));
    sum_w += w;

    z.coeffRef(i - 1) += sum_w;
    z.coeffRef(i) -= w * n;
  }

  return z;
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined using a modified version of the
 * the inverse of the isometric log ratio transform (ILR).
 * See:
 * Egozcue, Juan Jose; Pawlowsky-Glahn, Vera; Mateu-Figueras, Gloria;
 * Barcelo-Vidal, Carles (2003), "Isometric logratio transformations for
 * compositional data analysis", Mathematical Geology, 35 (3): 279–300,
 * doi:10.1023/A:1023818214614, S2CID 122844634
 *
 * This implementation is closer to the description of the same using "pivot
 * coordinates" in
 * Filzmoser, P., Hron, K., Templ, M. (2018). Geometrical Properties of
 * Compositional Data. In: Applied Compositional Data Analysis. Springer Series
 * in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-319-96422-5_3
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @param lp unused
 * @return Zero-sum vector of dimensionality K.
 */
template <typename Vec, require_eigen_col_vector_t<Vec>* = nullptr,
          require_not_st_var<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_constrain(const Vec& y,
                                               value_type_t<Vec>& lp) {
  return sum_to_zero_constrain(y);
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined using a modified version of
 * the inverse of the isometric log ratio transform (ILR).
 * See:
 * Egozcue, Juan Jose; Pawlowsky-Glahn, Vera; Mateu-Figueras, Gloria;
 * Barcelo-Vidal, Carles (2003), "Isometric logratio transformations for
 * compositional data analysis", Mathematical Geology, 35 (3): 279–300,
 * doi:10.1023/A:1023818214614, S2CID 122844634
 *
 * This implementation is closer to the description of the same using "pivot
 * coordinates" in
 * Filzmoser, P., Hron, K., Templ, M. (2018). Geometrical Properties of
 * Compositional Data. In: Applied Compositional Data Analysis. Springer Series
 * in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-319-96422-5_3
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Jacobian unused
 * @tparam Vec A type inheriting from `Eigen::DenseBase` or a `var_value` with
 *  inner type inheriting from `Eigen::DenseBase` with compile time dynamic rows
 *  and 1 column
 * @param[in] y free vector
 * @param[in, out] lp unused
 * @return Zero-sum vector of dimensionality one greater than `y`
 */
template <bool Jacobian, typename Vec, require_not_std_vector_t<Vec>* = nullptr>
inline plain_type_t<Vec> sum_to_zero_constrain(const Vec& y,
                                               return_type_t<Vec>& lp) {
  return sum_to_zero_constrain(y);
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined using a modified version of
 * the inverse of the isometric log ratio transform (ILR).
 * See:
 * Egozcue, Juan Jose; Pawlowsky-Glahn, Vera; Mateu-Figueras, Gloria;
 * Barcelo-Vidal, Carles (2003), "Isometric logratio transformations for
 * compositional data analysis", Mathematical Geology, 35 (3): 279–300,
 * doi:10.1023/A:1023818214614, S2CID 122844634
 *
 * This implementation is closer to the description of the same using "pivot
 * coordinates" in
 * Filzmoser, P., Hron, K., Templ, M. (2018). Geometrical Properties of
 * Compositional Data. In: Applied Compositional Data Analysis. Springer Series
 * in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-319-96422-5_3
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Jacobian unused
 * @tparam Vec A standard vector with inner type inheriting from
 * `Eigen::DenseBase` or a `var_value` with inner type inheriting from
 * `Eigen::DenseBase` with compile time dynamic rows and 1 column
 * @param[in] y free vector
 * @param[in, out] lp unused
 * @return Zero-sum vectors of dimensionality one greater than `y`
 */
template <bool Jacobian, typename T, require_std_vector_t<T>* = nullptr>
inline auto sum_to_zero_constrain(const T& y, return_type_t<T>& lp) {
  return apply_vector_unary<T>::apply(
      y, [](auto&& v) { return sum_to_zero_constrain(v); });
}

}  // namespace math
}  // namespace stan

#endif
