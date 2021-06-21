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
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<is_fvar, Ta, Tb>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, Ta> grad_a(a.size());
  promote_scalar_t<fvar_inner_t, Tb> grad_b(b.size());
  fvar_inner_t grad_z;
  grad_pFq(grad_a, grad_b, grad_z, a.val(), b.val(), z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(a.val(), b.val(), z.val()),
    dot_product(a.d(), grad_a) + dot_product(b.d(), grad_b) + z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_f, vector_f, double)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<is_fvar, Ta, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, Ta> grad_a(a.size());
  promote_scalar_t<fvar_inner_t, Tb> grad_b(b.size());
  grad_pFq_pq(grad_a, grad_b, a.val(), b.val(), z);
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(a.val(), b.val(), z),
    dot_product(a.d(), grad_a) + dot_product(b.d(), grad_b));
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_f, vector_d, double)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_eigen_vector_vt<is_fvar, Ta>* = nullptr,
          require_eigen_vector_vt<std::is_arithmetic, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, plain_type_t<Ta>> grad_a(a.size());
  grad_pFq_p(grad_a, a.val(), b, z);
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(a.val(), b, z),
    dot_product(a.d(), grad_a));
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_f, vector_d, fvar<T>)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<is_fvar, Ta>* = nullptr,
          require_all_eigen_vector_vt<std::is_arithmetic, Tb>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, Ta> grad_a(a.size());
  fvar_inner_t grad_z;
  grad_pFq_pz(grad_a, grad_z, a.val(), b, z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(a.val(), b, z.val()),
    dot_product(a.d(), grad_a) + z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_d, fvar<T>)
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Ta, Tb>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  fvar_inner_t grad_z;
  grad_pFq_z(grad_z, a, b, z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(a, b, z.val()),
    z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_f, fvar<T>)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Ta>* = nullptr,
          require_all_eigen_vector_vt<is_fvar, Tb>* = nullptr,
          require_fvar_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, plain_type_t<Tb>> grad_b(b.size());
  fvar_inner_t grad_z;
  grad_pFq_qz(grad_b, grad_z, a, b.val(), z.val());
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(a, b.val(), z.val()),
    dot_product(b.d(), grad_b) + z.d_ * grad_z);
}

/**
 * Returns the generalised hypergeometric function applied to the
 * input arguments.
 * 
 * Specialization for hypergeometric_pFq(vector_d, vector_f, double)
 * where vector_f = Eigen::Matrix<fvar<T>, -1, 1>
 *
 * @tparam Ta Type of a argument
 * @tparam Tb Type of b argument
 * @tparam Tz Type of z argument
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @return Generalised hypergeometric function
 */
template <typename Ta, typename Tb, typename Tz,
          require_all_eigen_vector_vt<std::is_arithmetic, Ta>* = nullptr,
          require_all_eigen_vector_vt<is_fvar, Tb>* = nullptr,
          require_arithmetic_t<Tz>* = nullptr>
inline return_type_t<Ta, Tb, Tz> hypergeometric_pFq(const Ta& a, const Tb& b,
                                                    const Tz& z) {
  using fvar_t = return_type_t<Ta, Tb, Tz>;
  using fvar_inner_t = typename fvar_t::Scalar;
  promote_scalar_t<fvar_inner_t, plain_type_t<Tb>> grad_b(b.size());
  grad_pFq_q(grad_b, a, b.val(), z);
  return fvar<fvar_inner_t>(
    hypergeometric_pFq(a, b.val(), z),
    dot_product(b.d(), grad_b));
}

}  // namespace math
}  // namespace stan
#endif
