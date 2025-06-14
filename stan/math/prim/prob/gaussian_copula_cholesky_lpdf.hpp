#ifndef STAN_MATH_PRIM_PROB_GAUSSIAN_COPULA_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_GAUSSIAN_COPULA_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inv_Phi.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/rep_vector.hpp>
#include <stan/math/prim/fun/vector_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_mvt.hpp>
#include <stan/math/prim/prob/multi_normal_cholesky_lpdf.hpp>
#include <stan/math/prim/prob/std_normal_lpdf.hpp>

namespace stan {
namespace math {
namespace internal {

template <typename Ty_elem, typename Ttuple, std::size_t... I>
auto apply_cdfs_elem_impl(Ty_elem&& y_elem, Ttuple&& cdf_tuple, std::index_sequence<I...>){
  auto&& cdf_func = std::get<0>(cdf_tuple);
  return cdf_func(
    std::forward<Ty_elem>(y_elem),
    std::get<I + 1>(std::forward<Ttuple>(cdf_tuple))...
  );
}


template <typename Ty, typename Ttuple, std::size_t... I>
auto apply_cdfs_impl(Ty&& y, Ttuple&& cdf_tuple, std::index_sequence<I...>){
  plain_type_t<Ty> res(y.size());

  // Use initializer-list expansion to assign the result of each CDF
  // to the respective element in the result vector, as we need a constant
  // expression for indexing the tuple of CDF functors
  static_cast<void>(std::initializer_list<int>{(
    res[I] = apply_cdfs_elem_impl(
      std::forward<decltype(y[I])>(y[I]),
      std::get<I>(cdf_tuple),
        std::make_index_sequence<
          // Using size - 1 as the first element is the functor to apply
          std::tuple_size<
            std::remove_reference_t<
              decltype(std::get<I>(cdf_tuple))>>{} - 1>{}),
    0)...});

  return res;
}

template <typename Ty, typename Ttuple>
auto apply_cdfs(Ty&& y, Ttuple&& cdf_tuple){
  return apply_cdfs_impl(
    std::forward<Ty>(y),
    std::forward<Ttuple>(cdf_tuple),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Ttuple>>{}>{}
  );
}
}

/*
  vectors ~ gaussian_copula_cholesky(cdf_funs_tuple, chol)
  (cdf_fun, ...)
*/

template <typename T_y, typename T_cdf_fun_tuple, typename T_chol,
          typename T_return = return_type_t<T_y, T_cdf_fun_tuple, T_chol>>
T_return gaussian_copula_cholesky_lpdf(
  const T_y& y, const T_cdf_fun_tuple& cdf_fun_tuple, const T_chol chol) {
  static constexpr const char* function = "gaussian_copula_cholesky_lpdf";

  using T_y_ref = ref_type_t<T_y>;
  using T_chochol_ref = ref_type_t<T_chol>;
  using T_cdf_ref = ref_type_t<T_cdf_fun_tuple>;

  check_consistent_sizes_mvt(function, "y", y, "cdf_fun_tuple", cdf_fun_tuple);
  const size_t size_mvt_y = size_mvt(y);
  const size_t size_mvt_cdf_tuple = size_mvt(cdf_fun_tuple);
  if (size_mvt_y == 0) {
    return 0;
  }
  T_y_ref y_ref = y;
  T_chochol_ref chol_ref = chol;
  T_cdf_ref cdf_tuple_ref = cdf_fun_tuple;

  vector_seq_view<T_y_ref> y_vec(y_ref);
  vector_seq_view<T_cdf_ref> cdf_tuple_vec(cdf_tuple_ref);
  const size_t size_vec = max_size_mvt(y, cdf_fun_tuple);

  const int size_y = math::size(y_vec[0]);
  const int size_cdf_tuple = math::size(cdf_tuple_vec[0]);

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
  for (size_t i = 1; i < size_mvt_cdf_tuple; i++) {
    check_size_match(function,
                     "Size of one of the vectors of "
                     "the CDF functor tuple",
                     math::size(cdf_tuple_vec[i]),
                     "Size of the first vector of the "
                     "location variable",
                     size_cdf_tuple);
  }

  check_size_match(function, "Size of random variable", size_y,
                   "size of CDF functor tuple", size_cdf_tuple);
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
    q[i] = to_vector(inv_Phi(internal::apply_cdfs(y_vec[i], cdf_tuple_vec[i])));
    lp -= std_normal_lpdf(q[i]);
  }
  const std::vector<Eigen::VectorXd> zero_vec(size_vec, rep_vector(0, size_y));
  lp += multi_normal_cholesky_lpdf(q, zero_vec, chol_ref);
  return lp;
}

}  // namespace math
}  // namespace stan
#endif
