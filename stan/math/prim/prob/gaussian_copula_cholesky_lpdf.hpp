#ifndef STAN_MATH_PRIM_PROB_GAUSSIAN_COPULA_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GAUSSIAN_COPULA_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/rep_vector.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp>
#include <stan/math/prim/prob/std_normal_lpdf.hpp>
#include <stan/math/prim/prob/std_normal_log_qf.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename Ty_elem, typename Ttuple, std::size_t... I>
auto apply_lcdfs_elem_impl(
  Ty_elem&& y_elem, Ttuple&& lcdf_tuple, std::index_sequence<I...>) {
  auto&& lcdf_func = std::get<0>(lcdf_tuple);
  return lcdf_func(
    std::forward<Ty_elem>(y_elem),
    std::get<I + 1>(std::forward<Ttuple>(lcdf_tuple))...
  );
}

template <typename Ty, typename Ttuple>
auto apply_lcdfs(Ty&& y, Ttuple&& lcdf_tuple) {
  return index_apply<std::tuple_size<std::remove_reference_t<Ttuple>>{}>(
    [&y, &lcdf_tuple](auto... Is) {
      return std::make_tuple(
        apply_lcdfs_elem_impl(
          y[Is], std::get<Is>(lcdf_tuple), std::make_index_sequence<
            std::tuple_size<std::remove_reference_t<decltype(std::get<Is>(lcdf_tuple))>>{} - 1>{}
        )...
      );
    });
}
}


/** \ingroup multivar_dists
 * The log of the Gaussian copula density for the given y and
 * a Cholesky factor L of the variance matrix.
 *
 * As the Gaussian copula requires inputs to be in the unit interval,
 * users are required to provide a tuple (of the same size as y) where each
 * element is a tuple containing a functor for calculating the LCDF and any
 * additional parameters required by the LCDF functor.
 *
 * This version of the function is vectorized on y and the tuple of LCDF
 * functor tuples.
 *
 * @param y A scalar vector
 * @param lcdf_fun_tuple A tuple of tuples, where each inner tuple follows the
 *   structure {functor, param1, param2, ...}, where the functor has signature
 *   (y[i], param1, param2, ...) and returns the log CDF for the given y[i].
 * @param chol The Cholesky decomposition of a variance matrix
 * of the multivariate normal distribution
 * @return The log of the gaussian copula density.
 * @throw std::domain_error if LL' is not square, not symmetric,
 * or not semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_lcdf_fun_tuple Type of tuple of LCDF functor tuples.
 * @tparam T_chol Type of cholesky factor.
 */
template <bool propto, typename T_y, typename T_lcdf_fun_tuple, typename T_chol>
auto gaussian_copula_cholesky_lpdf(
  const T_y& y, const T_lcdf_fun_tuple& lcdf_fun_tuple, const T_chol chol) {
  static constexpr const char* function = "gaussian_copula_cholesky_lpdf";

  using T_y_ref = ref_type_t<T_y>;
  using T_chol_ref = ref_type_t<T_chol>;
  using T_lcdf_ref = ref_type_t<T_lcdf_fun_tuple>;

  check_consistent_sizes_mvt(function, "y", y, "lcdf_fun_tuple", lcdf_fun_tuple);
  const size_t size_mvt_y = size_mvt(y);
  const size_t size_mvt_lcdf_tuple = size_mvt(lcdf_fun_tuple);
  T_y_ref y_ref = y;
  T_chol_ref chol_ref = chol;
  T_lcdf_ref lcdf_tuple_ref = lcdf_fun_tuple;

  vector_seq_view<T_y_ref> y_vec(y_ref);
  vector_seq_view<T_lcdf_ref> lcdf_tuple_vec(lcdf_tuple_ref);
  using T_return = return_type_t<T_y, decltype(internal::apply_lcdfs(y_vec[0], lcdf_tuple_vec[0])), T_chol>;

  if (size_mvt_y == 0) {
    return T_return(0);
  }
  const size_t size_vec = max_size_mvt(y, lcdf_fun_tuple);

  const int size_y = math::size(y_vec[0]);
  const int size_lcdf_tuple = math::size(lcdf_tuple_vec[0]);

  // check size consistency of all random variables y
  for (size_t i = 1; i < size_mvt_y; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the random variable",
                     math::size(y_vec[i]),
                     "Size of the first vector of the "
                     "random variable",
                     size_y);
  }
  // check size consistency of all CDF tuples
  for (size_t i = 1; i < size_mvt_lcdf_tuple; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the LCDF functor tuple",
                     math::size(lcdf_tuple_vec[i]),
                     "Size of the first vector of the "
                     "location variable",
                     size_lcdf_tuple);
  }

  check_size_match(function, "Size of random variable", size_y,
                   "size of LCDF functor tuple", size_lcdf_tuple);
  check_size_match(function, "Size of random variable", size_y,
                   "rows of covariance parameter", chol.rows());
  check_size_match(function, "Size of random variable", size_y,
                   "columns of covariance parameter", chol.cols());

  for (size_t i = 0; i < size_vec; i++) {
    check_not_nan(function, "Random variable", y_vec[i]);
  }
  check_cholesky_factor(function, "Cholesky decomposition of a variance matrix",
                        chol_ref);

  promote_scalar_t<T_return, std::vector<Eigen::VectorXd>> q(size_vec);
  T_return lp(0);
  for (size_t i = 0; i < size_vec; i++) {
    const auto& y_i = y_vec[i];
    const auto& func_i = lcdf_tuple_vec[i];

    const auto& res = internal::apply_lcdfs(y_i, func_i);
    const auto& u = index_apply<std::tuple_size<std::remove_reference_t<decltype(res)>>{}>(
      [&res, size_y](auto... Is) {
        Eigen::Matrix<T_return, Eigen::Dynamic, 1> u_inner(size_y);
        static_cast<void>(std::initializer_list<int>{
          (u_inner[Is] = std::get<Is>(res), 0)...
        });
        return u_inner;
      });
    check_bounded(function, "LCDF-transformed inputs", u, NEGATIVE_INFTY, 0);
    q[i] = std_normal_log_qf(u);
    lp -= std_normal_lpdf<propto>(q[i]);
  }
  const std::vector<Eigen::VectorXd> zero_vec(size_vec, rep_vector(0, size_y));
  lp += multi_normal_cholesky_lpdf<propto>(q, zero_vec, chol_ref);
  return lp;
}

template <typename T_y, typename T_lcdf_fun_tuple, typename T_chol>
auto gaussian_copula_cholesky_lpdf(
  const T_y& y, const T_lcdf_fun_tuple& lcdf_fun_tuple, const T_chol chol) {
  return gaussian_copula_cholesky_lpdf<false>(y, lcdf_fun_tuple, chol);
}
}  // namespace math
}  // namespace stan
#endif
