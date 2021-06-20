#ifndef STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP
#define STAN_MATH_FWD_FUN_HYPERGEOMETRIC_PFQ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/hypergeometric_pFq.hpp>
#include <stan/math/prim/fun/grad_pFq.hpp>

namespace stan {
namespace math {

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_f, vector_f, fvar<T>)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Tp Type of p argument
 * @tparam Tq Type of q argument
 * @tparam Tz Type of z argument
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<is_fvar, Tp, Tq>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Tp, Tq, Tz> hypergeometric_pFq(const Tp& p, const Tq& q,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Tp, Tq, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, Tp> grad_p(p.size());
  promote_scalar_t<fvar_inner_t, Tq> grad_q(q.size());
  fvar_inner_t grad_z;
  grad_pFq(grad_p, grad_q, grad_z, p.val(), q.val(), z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(p.val(), q.val(), z.val()),
    dot_product(p.d(), grad_p) + dot_product(q.d(), grad_q) + z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_f, vector_f, double)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Tp Type of p argument
 * @tparam Tq Type of q argument
 * @tparam Tz Type of z argument
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<is_fvar, Tp, Tq>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline return_type_t<Tp, Tq, Tz> hypergeometric_pFq(const Tp& p, const Tq& q,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Tp, Tq, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, Tp> grad_p(p.size());
  promote_scalar_t<fvar_inner_t, Tq> grad_q(q.size());
  grad_pFq_pq(grad_p, grad_q, p.val(), q.val(), z);
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(p.val(), q.val(), z),
    dot_product(p.d(), grad_p) + dot_product(q.d(), grad_q));
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_f, vector_d, double)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Tp Type of p argument
 * @tparam Tq Type of q argument
 * @tparam Tz Type of z argument
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Tp, typename Tq, typename Tz,
          require_eigen_vector_vt<is_fvar, Tp>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, Tq>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline return_type_t<Tp, Tq, Tz> hypergeometric_pFq(const Tp& p, const Tq& q,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Tp, Tq, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, plain_type_t<Tp>> grad_p(p.size());
  grad_pFq_p(grad_p, p.val(), q, z);
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(p.val(), q, z),
    dot_product(p.d(), grad_p));
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_f, vector_d, fvar<T>)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Tp Type of p argument
 * @tparam Tq Type of q argument
 * @tparam Tz Type of z argument
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<is_fvar, Tp>* = nullptr,
          require_all_eigen_vector_vt<std::is_arithmetic, Tq>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Tp, Tq, Tz> hypergeometric_pFq(const Tp& p, const Tq& q,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Tp, Tq, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, Tp> grad_p(p.size());
  fvar_inner_t grad_z;
  grad_pFq_pz(grad_p, grad_z, p.val(), q, z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(p.val(), q, z.val()),
    dot_product(p.d(), grad_p) + z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_d, fvar<T>)
 *
 * @tparam Tp Type of p argument
 * @tparam Tq Type of q argument
 * @tparam Tz Type of z argument
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Tp, Tq>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Tp, Tq, Tz> hypergeometric_pFq(const Tp& p, const Tq& q,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Tp, Tq, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  fvar_inner_t grad_z;
  grad_pFq_z(grad_z, p, q, z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(p, q, z.val()),
    z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_f, fvar<T>)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Tp Type of p argument
 * @tparam Tq Type of q argument
 * @tparam Tz Type of z argument
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Tp>* = nullptr,
          require_all_eigen_vector_vt<is_fvar, Tq>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Tp, Tq, Tz> hypergeometric_pFq(const Tp& p, const Tq& q,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Tp, Tq, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, plain_type_t<Tq>> grad_q(q.size());
  fvar_inner_t grad_z;
  grad_pFq_qz(grad_q, grad_z, p, q.val(), z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(p, q.val(), z.val()),
    dot_product(q.d(), grad_q) + z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_f, double)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Tp Type of p argument
 * @tparam Tq Type of q argument
 * @tparam Tz Type of z argument
 * @param[in] p p
 * @param[in] q q
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Tp, typename Tq, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Tp>* = nullptr,
          require_all_eigen_vector_vt<is_fvar, Tq>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline return_type_t<Tp, Tq, Tz> hypergeometric_pFq(const Tp& p, const Tq& q,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Tp, Tq, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, plain_type_t<Tq>> grad_q(q.size());
  grad_pFq_q(grad_q, p, q.val(), z);
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(p, q.val(), z),
    dot_product(q.d(), grad_q));
}

}  // namespace math
}  // namespace stan
#endif