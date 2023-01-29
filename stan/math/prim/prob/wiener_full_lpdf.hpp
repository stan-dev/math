// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef STAN_MATH_PRIM_PROB_WIENER_FULL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WIENER_FULL_LPDF_HPP

#include <stan/math/prim/fun.hpp>

namespace stan {
namespace math {

namespace internal {

template <bool propto, typename T_y, typename T_alpha, typename T_delta,
          typename T_beta, typename T_t0, typename T_sv, typename T_sw,
          typename T_st0>
inline return_type_t<T_y, T_alpha, T_delta, T_beta, T_t0, T_sv, T_sw, T_st0>
wiener_full_prec_impl_lpdf(const char* function_name, const T_y& y,
                           const T_alpha& a, const T_delta& v, const T_beta& w,
                           const T_t0& t0, const T_sv& sv, const T_sw& sw,
                           const T_st0& st0, const double& prec) {
  using T_partials_return = partials_return_t<T_y, T_alpha, T_delta, T_beta,
                                              T_t0, T_sv, T_sw, T_st0>;
  using T_y_ref = ref_type_t<T_y>;
  using T_alpha_ref = ref_type_t<T_alpha>;
  using T_delta_ref = ref_type_t<T_delta>;
  using T_beta_ref = ref_type_t<T_beta>;
  using T_t0_ref = ref_type_t<T_t0>;
  using T_sv_ref = ref_type_t<T_sv>;
  using T_sw_ref = ref_type_t<T_sw>;
  using T_st0_ref = ref_type_t<T_st0>;

  static constexpr double LOG_FOUR = LOG_TWO + LOG_TWO;
  static constexpr double LOG_POINT1 = -1;

  check_consistent_sizes(function_name, "Random variable", y,
                         "Boundary separation", a, "Drift rate", v,
                         "A-priori bias", w, "Nondecision time", t0,
                         "Inter-trial variability in drift rate", sv,
                         "Inter-trial variability in A-priori bias", sw,
                         "Inter-trial variability in Nondecision time", st0);

  check_consistent_size(function_name, "Random variable", y, 1);
  check_consistent_size(function_name, "Boundary separation", a, 1);
  check_consistent_size(function_name, "Drift rate", v, 1);
  check_consistent_size(function_name, "A-priori bias", w, 1);
  check_consistent_size(function_name, "Nondecision time", t0, 1);
  check_consistent_size(function_name, "Inter-trial variability in drift rate",
                        sv, 1);
  check_consistent_size(function_name,
                        "Inter-trial variability in A-priori bias", sw, 1);
  check_consistent_size(function_name,
                        "Inter-trial variability in Nondecision time", st0, 1);

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = a;
  T_delta_ref delta_ref = v;
  T_beta_ref beta_ref = w;
  T_t0_ref t0_ref = t0;
  T_sv_ref sv_ref = sv;
  T_sw_ref sw_ref = sw;
  T_st0_ref st0_ref = st0;

  check_positive_finite(function_name, "Random variable", value_of(y_ref));
  check_positive_finite(function_name, "Boundary separation",
                        value_of(alpha_ref));
  check_finite(function_name, "Drift rate", value_of(delta_ref));
  check_less(function_name, "A-priori bias", value_of(beta_ref), 1);
  check_greater(function_name, "A-priori bias", value_of(beta_ref), 0);
  check_nonnegative(function_name, "Nondecision time", value_of(t0_ref));
  check_finite(function_name, "Nondecision time", value_of(t0_ref));
  check_nonnegative(function_name, "Inter-trial variability in drift rate",
                    value_of(sv_ref));
  check_finite(function_name, "Inter-trial variability in drift rate",
               value_of(sv_ref));
  check_bounded(function_name, "Inter-trial variability in A-priori bias",
                value_of(sw_ref), 0, 1);
  check_nonnegative(function_name,
                    "Inter-trial variability in Nondecision time",
                    value_of(st0_ref));
  check_finite(function_name, "Inter-trial variability in Nondecision time",
               value_of(st0_ref));

  if (size_zero(y, a, v, w, t0) || size_zero(sv, sw, st0))
    return 0;

  size_t N = max_size(y, a, v, w, t0, sv, sw, st0);
  if (!N)
    return 0;

  scalar_seq_view<T_y_ref> y_vec(y_ref);
  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_delta_ref> delta_vec(delta_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  scalar_seq_view<T_t0_ref> t0_vec(t0_ref);
  scalar_seq_view<T_sv_ref> sv_vec(sv_ref);
  scalar_seq_view<T_sw_ref> sw_vec(sw_ref);
  scalar_seq_view<T_st0_ref> st0_vec(st0_ref);
  size_t N_y_t0 = max_size(y, t0, st0);

  for (size_t i = 0; i < N_y_t0; ++i) {
    if (y_vec[i] <= t0_vec[i]) {
      std::stringstream msg;
      msg << ", but must be greater than nondecision time = " << t0_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name, "Random variable", y_vec[i], " = ",
                         msg_str.c_str());
    }
  }
  size_t N_beta_sw = max_size(w, sw);
  for (size_t i = 0; i < N_beta_sw; ++i) {
    if (beta_vec[i] - .5 * sw_vec[i] <= 0) {
      std::stringstream msg;
      msg << ", but must be smaller than 2*(A-priori bias) = "
          << 2 * beta_vec[i];
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
    if (beta_vec[i] + .5 * sw_vec[i] >= 1) {
      std::stringstream msg;
      msg << ", but must be smaller than 2*(1-A-priori bias) = "
          << 2 * (1 - beta_vec[i]);
      std::string msg_str(msg.str());
      throw_domain_error(function_name,
                         "Inter-trial variability in A-priori bias", sw_vec[i],
                         " = ", msg_str.c_str());
    }
  }

  if (!include_summand<propto, T_y, T_alpha, T_delta, T_beta, T_t0, T_sv, T_sw,
                       T_st0>::value)
    return 0;

  // abstol_wiener5 z.B. 1e-12 on normal scale
  // T_partials_return error_bound = 1e-6;
  // alpha in Skizze
  T_partials_return error_bound_dens = 1e-6;  // precision for
  // density (fixed here, test with 1e-6 and 1e-2)
  T_partials_return lerror_bound_dens = log(error_bound_dens);
  T_partials_return error_bound = prec;  // precision for
  // derivatives (controllable by user)
  T_partials_return lerror_bound = log(error_bound);  // log(alpha)
  T_partials_return abstol = 0.0;
  T_partials_return reltol = .9 * error_bound;  // eps_rel(Integration)
  T_partials_return abstol_wiener5 = 1e-12;     // eps_abs(wiener5)
  T_partials_return labstol_wiener5 = log(abstol_wiener5);
  // log(eps_abs(wiener5)
  int Meval = 6000;
  T_partials_return dens = 0.0;
  T_partials_return ld = 0.0;
  operands_and_partials<T_y_ref, T_alpha_ref, T_delta_ref, T_beta_ref, T_t0_ref,
                        T_sv_ref, T_sw_ref, T_st0_ref>
      ops_partials(y_ref, alpha_ref, delta_ref, beta_ref, t0_ref, sv_ref,
                   sw_ref, st0_ref);

  // calculate density and partials
  for (size_t i = 0; i < N; i++) {
    // Calculate 4-parameter model without inter-trial variabilities (if
    // sv_vec[i] == 0) or 5-parameter model with inter-trial variability in
    // drift rate (if sv_vec[i] != 0)
    if (sw_vec[i] == 0 && st0_vec[i] == 0) {
      const T_partials_return y_val = y_vec.val(i);
      const T_partials_return alpha_val = alpha_vec.val(i);
      const T_partials_return delta_val = delta_vec.val(i);
      const T_partials_return beta_val = beta_vec.val(i);
      const T_partials_return t0_val = t0_vec.val(i);
      const T_partials_return sv_val = sv_vec.val(i);

      dens = internal::dwiener5(y_val - t0_val, alpha_val, delta_val, beta_val,
                                sv_val, labstol_wiener5);
      if (labstol_wiener5 > fabs(dens) + lerror_bound_dens - LOG_TWO)
        dens = internal::dwiener5(y_val - t0_val, alpha_val, delta_val,
                                  beta_val, sv_val,
                                  fabs(dens) + lerror_bound_dens - LOG_TWO);
      ld += dens;

      // computation of derivative for t and precision check in order to give
      // the value as deriv_t to edge1 and as -deriv_t to edge5
      T_partials_return deriv_t
          = internal::dtdwiener5(y_val - t0_val, alpha_val, delta_val, beta_val,
                                 sv_val, labstol_wiener5);
      if (labstol_wiener5 > log(fabs(deriv_t)) + dens + lerror_bound - LOG_TWO)
        deriv_t = internal::dtdwiener5(
            y_val - t0_val, alpha_val, delta_val, beta_val, sv_val,
            log(fabs(deriv_t)) + dens + lerror_bound - LOG_FOUR);

      // computation of derivatives and precision checks
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[i] = deriv_t;
      }
      if (!is_constant_all<T_alpha>::value) {
        T_partials_return deriv_a
            = internal::dadwiener5(y_val - t0_val, alpha_val, delta_val,
                                   beta_val, sv_val, labstol_wiener5, 0);
        if (labstol_wiener5
            > log(fabs(deriv_a)) + dens + lerror_bound - LOG_TWO)
          deriv_a = internal::dadwiener5(
              y_val - t0_val, alpha_val, delta_val, beta_val, sv_val,
              log(fabs(deriv_a)) + dens + lerror_bound - LOG_FOUR, 0);
        ops_partials.edge2_.partials_[i] = deriv_a;
      }
      if (!is_constant_all<T_delta>::value) {
        ops_partials.edge3_.partials_[i] = internal::dvdwiener5(
            y_val - t0_val, alpha_val, delta_val, beta_val, sv_val);
      }
      if (!is_constant_all<T_beta>::value) {
        T_partials_return deriv_w
            = internal::dwdwiener5(y_val - t0_val, alpha_val, delta_val,
                                   beta_val, sv_val, labstol_wiener5, 0);
        if (labstol_wiener5
            > log(fabs(deriv_w)) + dens + lerror_bound - LOG_TWO)
          deriv_w = internal::dwdwiener5(
              y_val - t0_val, alpha_val, delta_val, beta_val, sv_val,
              log(fabs(deriv_w)) + dens + lerror_bound - LOG_FOUR, 0);
        ops_partials.edge4_.partials_[i] = deriv_w;
      }
      if (!is_constant_all<T_t0>::value) {
        ops_partials.edge5_.partials_[i] = -deriv_t;
      }
      if (!is_constant_all<T_sv>::value) {
        ops_partials.edge6_.partials_[i] = internal::dsvdwiener5(
            y_val - t0_val, alpha_val, delta_val, beta_val, sv_val);
      }
      if (!is_constant_all<T_sw>::value) {
        ops_partials.edge7_.partials_[i] = 0;
      }
      if (!is_constant_all<T_st0>::value) {
        ops_partials.edge8_.partials_[i] = 0;
      }
      // Calculate 6-, or 7-parameter model
    } else {
      const T_partials_return y_val = y_vec.val(i);
      const T_partials_return alpha_val = alpha_vec.val(i);
      const T_partials_return delta_val = delta_vec.val(i);
      const T_partials_return beta_val = beta_vec.val(i);
      const T_partials_return t0_val = t0_vec.val(i);
      const T_partials_return sv_val = sv_vec.val(i);
      const T_partials_return sw_val = sw_vec.val(i);
      const T_partials_return st0_val = st0_vec.val(i);
      internal::my_params<
          T_partials_return, T_partials_return, T_partials_return,
          T_partials_return, T_partials_return, T_partials_return,
          T_partials_return, T_partials_return, T_partials_return>
          params = {y_val,    alpha_val, delta_val,
                    beta_val, t0_val,    sv_val,
                    sw_val,   st0_val,   labstol_wiener5 - LOG_TWO};
      int dim = (sw_val != 0) + (st0_val != 0);
      check_positive(function_name,
                     "(Inter-trial variability in A-priori bias) + "
                     "(Inter-trial variability in nondecision time)",
                     dim);

      std::vector<T_partials_return> xmin(dim, 0);
      std::vector<T_partials_return> xmax(dim, 1);

      if (st0_val)
        xmax[dim - 1] = fmin(1.0, (y_val - t0_val) / st0_val);

      dens = hcubature(internal::int_ddiff<std::vector<double>, void*>, &params,
                       dim, xmin, xmax, Meval, abstol, reltol / 2);
      if (labstol_wiener5
          > fabs(log(dens)) + LOG_POINT1 + lerror_bound_dens - LOG_TWO) {
        T_partials_return new_error
            = LOG_POINT1 + lerror_bound_dens - LOG_TWO + log(fabs(dens));
        internal::my_params<
            T_partials_return, T_partials_return, T_partials_return,
            T_partials_return, T_partials_return, T_partials_return,
            T_partials_return, T_partials_return, T_partials_return>
            params_new_error = {y_val,  alpha_val, delta_val, beta_val, t0_val,
                                sv_val, sw_val,    st0_val,   new_error};
        dens = hcubature(internal::int_ddiff<std::vector<double>, void*>,
                         &params_new_error, dim, xmin, xmax, Meval, abstol,
                         reltol);
      }
      ld += log(dens);

      // computation of derivative for t and precision check in order to give
      // the value as deriv_t to edge1 and as -deriv_t to edge5
      T_partials_return deriv_t_7;
      deriv_t_7
          = 1 / dens
            * hcubature(internal::int_dtddiff<std::vector<double>, void*>,
                        &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
      if (labstol_wiener5
          > log(fabs(deriv_t_7)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
        T_partials_return new_error
            = LOG_POINT1 + lerror_bound - LOG_FOUR + log(fabs(deriv_t_7));
        internal::my_params<
            T_partials_return, T_partials_return, T_partials_return,
            T_partials_return, T_partials_return, T_partials_return,
            T_partials_return, T_partials_return, T_partials_return>
            params_new_error = {y_val,  alpha_val, delta_val, beta_val, t0_val,
                                sv_val, sw_val,    st0_val,   new_error};
        deriv_t_7
            = 1 / dens
              * hcubature(internal::int_dtddiff<std::vector<double>, void*>,
                          &params_new_error, dim, xmin, xmax, Meval, abstol,
                          reltol / 2);
      }

      // computation of derivatives and precision checks
      T_partials_return deriv;
      if (!is_constant_all<T_y>::value) {
        ops_partials.edge1_.partials_[i] = deriv_t_7;
      }
      if (!is_constant_all<T_alpha>::value) {
        deriv
            = 1 / dens
              * hcubature(internal::int_daddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
        if (labstol_wiener5
            > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
          T_partials_return new_error
              = LOG_POINT1 + lerror_bound - LOG_FOUR + log(fabs(deriv));
          internal::my_params<
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return>
              params_new_error
              = {y_val,  alpha_val, delta_val, beta_val, t0_val,
                 sv_val, sw_val,    st0_val,   new_error};
          deriv = 1 / dens
                  * hcubature(internal::int_daddiff<std::vector<double>, void*>,
                              &params_new_error, dim, xmin, xmax, Meval, abstol,
                              reltol / 2);
        }
        ops_partials.edge2_.partials_[i] = deriv;
      }
      if (!is_constant_all<T_delta>::value) {
        deriv
            = 1 / dens
              * hcubature(internal::int_dvddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
        if (labstol_wiener5
            > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
          T_partials_return new_error
              = LOG_POINT1 + lerror_bound - LOG_TWO + log(fabs(deriv));
          internal::my_params<
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return>
              params_new_error
              = {y_val,  alpha_val, delta_val, beta_val, t0_val,
                 sv_val, sw_val,    st0_val,   new_error};
          deriv = 1 / dens
                  * hcubature(internal::int_dvddiff<std::vector<double>, void*>,
                              &params_new_error, dim, xmin, xmax, Meval, abstol,
                              reltol / 2);
        }
        ops_partials.edge3_.partials_[i] = deriv;
      }
      if (!is_constant_all<T_beta>::value) {
        deriv
            = 1 / dens
              * hcubature(internal::int_dwddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
        if (labstol_wiener5
            > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
          T_partials_return new_error
              = LOG_POINT1 + lerror_bound - LOG_FOUR + log(fabs(deriv));
          internal::my_params<
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return>
              params_new_error
              = {y_val,  alpha_val, delta_val, beta_val, t0_val,
                 sv_val, sw_val,    st0_val,   new_error};
          deriv = 1 / dens
                  * hcubature(internal::int_dwddiff<std::vector<double>, void*>,
                              &params_new_error, dim, xmin, xmax, Meval, abstol,
                              reltol / 2);
        }
        ops_partials.edge4_.partials_[i] = deriv;
      }
      if (!is_constant_all<T_t0>::value) {
        ops_partials.edge5_.partials_[i] = -deriv_t_7;
      }
      if (!is_constant_all<T_sv>::value) {
        deriv
            = 1 / dens
              * hcubature(internal::int_dsvddiff<std::vector<double>, void*>,
                          &params, dim, xmin, xmax, Meval, abstol, reltol / 2);
        if (labstol_wiener5
            > log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO) {
          T_partials_return new_error
              = LOG_POINT1 + lerror_bound - LOG_TWO + log(fabs(deriv));
          internal::my_params<
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return,
              T_partials_return, T_partials_return, T_partials_return>
              params_new_error
              = {y_val,  alpha_val, delta_val, beta_val, t0_val,
                 sv_val, sw_val,    st0_val,   new_error};
          deriv
              = 1 / dens
                * hcubature(internal::int_dsvddiff<std::vector<double>, void*>,
                            &params_new_error, dim, xmin, xmax, Meval, abstol,
                            reltol / 2);
        }
        ops_partials.edge6_.partials_[i] = deriv;
      }
      if (!is_constant_all<T_sw>::value) {
        if (sw_val == 0) {
          ops_partials.edge7_.partials_[i] = 0;
        } else {
          T_partials_return lower, upper, width, fl, fu;

          lower = beta_val - sw_val / 2;
          lower = (0 > lower) ? 0 : lower;
          upper = beta_val + sw_val / 2;
          upper = (1 < upper) ? 1 : upper;
          width = upper - lower;

          int dim_ = (st0_val != 0);
          if (dim_ == 0) {
            fl = exp(internal::dwiener5(y_val - t0_val, alpha_val, delta_val,
                                        lower, sv_val, labstol_wiener5));
            fu = exp(internal::dwiener5(y_val - t0_val, alpha_val, delta_val,
                                        upper, sv_val, labstol_wiener5));
            if (labstol_wiener5 > log(fabs(fl + fu) + log(sw_val) - LOG_TWO
                                      + lerror_bound - LOG_TWO)) {
              fl = exp(internal::dwiener5(y_val - t0_val, alpha_val, delta_val,
                                          lower, sv_val,
                                          lerror_bound - LOG_TWO + log(sw_val)
                                              - LOG_TWO + log(fabs(fl + fu))));
              fu = exp(internal::dwiener5(y_val - t0_val, alpha_val, delta_val,
                                          upper, sv_val,
                                          lerror_bound - LOG_TWO + log(sw_val)
                                              - LOG_TWO + log(fabs(fl + fu))));
            }
            deriv = 1 / width * 0.5 * (fl + fu);
          } else {
            internal::my_params2<
                T_partials_return, T_partials_return, T_partials_return,
                T_partials_return, T_partials_return, T_partials_return,
                T_partials_return, T_partials_return, T_partials_return,
                T_partials_return, T_partials_return, T_partials_return>
                params_sw
                = {y_val, alpha_val, delta_val, beta_val,
                   lower, upper,     t0_val,    sv_val,
                   0,     sw_val,    st0_val,   labstol_wiener5 - LOG_TWO};
            deriv = hcubature(
                internal::int_dswddiff<std::vector<double>, void*>, &params_sw,
                dim_, xmin, xmax, Meval, abstol, reltol / 2);
            if (labstol_wiener5
                > log(fabs(deriv) + LOG_POINT1 + lerror_bound - LOG_TWO)) {
              T_partials_return new_error
                  = log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO;
              internal::my_params2<
                  T_partials_return, T_partials_return, T_partials_return,
                  T_partials_return, T_partials_return, T_partials_return,
                  T_partials_return, T_partials_return, T_partials_return,
                  T_partials_return, T_partials_return, T_partials_return>
                  params_new_error_sw
                  = {y_val, alpha_val, delta_val, beta_val,
                     lower, upper,     t0_val,    sv_val,
                     0,     sw_val,    st0_val,   new_error};
              deriv = hcubature(
                  internal::int_dswddiff<std::vector<double>, void*>,
                  &params_new_error_sw, dim_, xmin, xmax, Meval, abstol,
                  reltol / 2);
            }
          }
          ops_partials.edge7_.partials_[i] = deriv / dens - 1 / sw_val;
        }
      }
      if (!is_constant_all<T_st0>::value) {
        T_partials_return f;
        if (st0_val == 0) {
          ops_partials.edge8_.partials_[i] = 0;
        } else if (y_val - (t0_val + st0_val) <= 0) {
          ops_partials.edge8_.partials_[i] = -1 / st0_val;
        } else {
          int dim_ = (sw_val != 0);
          if (dim_ == 0) {
            f = exp(internal::dwiener5(y_val - (t0_val + st0_val), alpha_val,
                                       delta_val, beta_val, sv_val,
                                       labstol_wiener5));
            if (labstol_wiener5 > log(fabs(f) + log(st0_val) - LOG_TWO
                                      + lerror_bound - LOG_TWO))
              f = exp(internal::dwiener5(y_val - (t0_val + st0_val), alpha_val,
                                         delta_val, beta_val, sv_val,
                                         lerror_bound - LOG_TWO + log(st0_val)
                                             - LOG_TWO + log(fabs(f))));
            deriv = 1 / st0_val * f;
          } else {
            T_partials_return new_error = labstol_wiener5 - LOG_TWO;
            internal::my_params3<
                T_partials_return, T_partials_return, T_partials_return,
                T_partials_return, T_partials_return, T_partials_return,
                T_partials_return, T_partials_return, T_partials_return,
                T_partials_return, T_partials_return>
                params_st = {y_val,    alpha_val, delta_val,
                             beta_val, t0_val,    t0_val + st0_val,
                             sv_val,   sw_val,    st0_val,
                             0,        new_error};
            deriv = hcubature(
                internal::int_dst0ddiff<std::vector<double>, void*>, &params_st,
                dim_, xmin, xmax, Meval, abstol, reltol / 2);
            if (labstol_wiener5
                > log(fabs(deriv) + LOG_POINT1 + lerror_bound - LOG_TWO)) {
              T_partials_return new_error
                  = log(fabs(deriv)) + LOG_POINT1 + lerror_bound - LOG_TWO;
              internal::my_params3<
                  T_partials_return, T_partials_return, T_partials_return,
                  T_partials_return, T_partials_return, T_partials_return,
                  T_partials_return, T_partials_return, T_partials_return,
                  T_partials_return, T_partials_return>
                  params_new_error_st = {y_val,    alpha_val, delta_val,
                                         beta_val, t0_val,    t0_val + st0_val,
                                         sv_val,   sw_val,    st0_val,
                                         0,        new_error};
              deriv = hcubature(
                  internal::int_dst0ddiff<std::vector<double>, void*>,
                  &params_new_error_st, dim_, xmin, xmax, Meval, abstol,
                  reltol / 2);
            }
          }
          ops_partials.edge8_.partials_[i] = -1 / st0_val + deriv / dens;
        }
      }
      std::vector<double>().swap(xmin);
      std::vector<double>().swap(xmax);
    }
  }

  return ops_partials.build(ld);
}
//-----------------------------------------------

}  // namespace internal

/** \ingroup prob_dists
 * The log of the first passage time density function for a (Wiener)
 * drift diffusion model with up to 7 parameters, where
 * \f$y\in \mathbb{R}_{+}\f$ is the reacion time, \f$a \in \mathbb{R}_{+}\f$
 * the boundary separation, \f$v \in \mathbb{R}\f$ the drifte rate,
 * \f$w \in (0, 1)\f$ the relative starting point (aka a-priori bias),
 * \f$t_0 \in \mathbb{R}_{\geq 0}\f$ the non-decision time, \f$s_v \in
 * \mathbb{R}_{\geq 0}\f$ the inter-trial variability of the drift rate,
 * \f$s_w \in [0, 1)\f$ the inter-trial  variability of the relative starting
 * point, and  \f$s_{t_0} \in \mathbb{R}_{\geq 0}\f$ the inter-trial variability
 * of
 * the non-decision time.
 *
 *
 \f{eqnarray*}{
 y &\sim& \text{wiener_full}(a,v,w,t_0,s_v,s_w,s_{t_0}) \\
 \log(p(y|a,v,w,t_0,s_v,s_w,s_{t_0})) &=& \log(\frac{1}{s_{t_0}}
 \int_{t_0}^{t_o + s_{t_0}} \frac{1}{s_{w}}\int_{w -0.5s_w}^{w + 0.5s_{w}}
 \int_{-\infty}^{\infty} p_3(y-\tau_0|a,\nu,\omega) \\
 &&\times\frac{1}{\sqrt{2\pi s_\nu^2}}
 \mathbb{e}^{-\frac{(\nu-v)^2}{2s_\nu^2}} \ d\nu \ d\omega \ d\tau_0) \\
 &=& \log(\frac{1}{s_{t_0}}
 \int_{t_0}^{t_o + s_{t_0}} \frac{1}{s_{w}}\int_{w -0.5s_w}^{w + 0.5s_{w}}
 M \times p_3(y-\tau_0|a,v,\omega) \ d\omega \ d\tau_0),
 \f}
 * where \f$M\f$ and \f$p_3()\f$ are defined, by using \f$t:=y-\tau_0\f$, as
 \f{eqnarray*}{
 M &:=& \frac{1}{\sqrt{1+s^2_v t}}
 \mathbb{e}^{av\omega+\frac{v^2t}{2}+\frac{s^2_v a^2
 \omega^2-2av\omega-v^2t}{2(1+s^2_vt)}} \text{ and} \\ p_3(t|a,v,w) &:=&
 \frac{1}{a^2} \mathbb{e}^{-a v w -\frac{v^2t}{2}} f(\frac{t}{a^2}|0,1,w), \f}
 * where  \f$f(t^*=\frac{t}{a^2}|0,1,w)\f$ has two forms
 \f{eqnarray*}{
 f_l(t^*|0,1,w) &=& \sum_{k=1}^{\infty} k\pi \mathbb{e}^{-\frac{k^2\pi^2t^*}{2}}
 \sin{(k \pi w)}\text{ and} \\
 f_s(t^*|0,1,w) &=& \sum_{k=-\infty}^{\infty} \frac{1}{\sqrt{2\pi (t^*)^3}}
 (w+2k) \mathbb{e}^{-\frac{(w+2k)^2}{2t^*}}, \f}
 * which are selected depending on the number of components \f$k\f$, needed to
 * guarantee a certain precision.
 *
 * Note that the parametrization for non-decision time and relative starting
 point is as follows:
 * \c t0 is the lower bound of the variability interval;
 * \c w is the mean of the variability interval.
 *
 * See \b Details below for more details on how to use \c wiener_full_lpdf().
 *
 * @tparam T_y type of scalar
 * @tparam T_alpha type of boundary
 * @tparam T_delta type of drift rate
 * @tparam T_beta type of relative starting point
 * @tparam T_t0 type of non-decision time
 * @tparam T_sv type of inter-trial variability of drift rate
 * @tparam T_sw type of inter-trial variability of relative starting point
 * @tparam T_st0 type of inter-trial variability of non-decision time
 *
 * @param y A scalar variable; the reaction time in seconds
 * @param a The boundary separation
 * @param v The drift rate
 * @param w The relative starting point
 * @param t0 The non-decision time
 * @param sv The inter-trial variability of the drift rate
 * @param sw The inter-trial variability of the relative starting point
 * @param st0 The inter-trial variability of the non-decision time
 * @return The log of the Wiener first passage time density with
 *  the specified arguments for upper boundary responses
 * @throw std::domain_error if non-decision time \c t0 is greater than reaction
 time \c y.
 * @throw std::domain_error if \c 1-sw/2 is smaller than or equal to \c w.
 * @throw std::domain_error if \c sw/2 is larger than or equal to \c w.
 *
 *
 *
 * **Details**
 *
 * The function can be called by
 @code
 target += wiener_full_lpdf(y, a, v, w, t0, sv, sw, st0);
 @endcode
 * or
 @code
 y ~ wiener_full(a, v, w, t0, sv, sw, st0);
 @endcode
 *
 * By default \c wiener_full_lpdf() gives the log of the
 * Wiener first-passage time probability density function for the \e upper
 response
 * boundary. To use the \e lower response boundary \c v and \c w must be changed
 to
 * \c -v and \c (1-w), respectively.
 *
 * \c sv, \c sw, \c st0, and \c t0 can be set to zero, indicating no inter-trial
 * variability in \f$v\f$, no inter-trial variability in \f$w\f$, no
 * inter-trial variability in \f$t_0\f$, and no non-decision time, respectively.
 * If \c t0 is zero, \c st0 must be zero as well. For example, when no
 inter-trial
 * variability for the relative starting point is needed one can write something
 like:
 @code
 target += wiener_full_lpdf(y, a, v, w, t0, sv, 0, st0)
 @endcode
 * If no inter-trial variability is needed at all one can write something like:
 @code
 target += wiener_full_lpdf(y, a, v, w, t0, 0, 0, 0)
 @endcode
 * If for some reason no non-decision time is assumed one can write something
 like:
 @code
 target += wiener_full_lpdf(y, a, v, w, 0, sv, sw, 0)
 @endcode
 *
 * To also control the precision in the estimation of the partial derivatives,
 * call the function \c wiener_full_prec_lpdf(), analogously:
 @code
 target += wiener_full__prec_lpdf(y, a, v, w, t0, sv, sw, st0, precision);
 @endcode
 *
 *
 * **References**
 * - Blurton, S. P., Kesselmeier, M., & Gondan, M. (2017). The first-passage
 * time distribution for the diffusion model with variable drift.
 * *Journal of Mathematical Psychology, 76*, 7–12.
 * https://doi.org/10.1016/j.jmp.2016.11.003
 * - Foster, K., & Singmann, H. (2021). Another Approximation of the
 First-Passage
 * Time Densities for the Ratcliff Diffusion Decision Model.
 * *arXiv preprint arXiv:2104.01902*
 * - Gondan, M., Blurton, S. P., & Kesselmeier, M. (2014). Even faster and even
 * more accurate first-passage time densities and distributions for the Wiener
 * diffusion model. *Journal of Mathematical Psychology, 60*, 20–22.
 * https://doi.org/10.1016/j.jmp.2014.05.002
 * - Hartmann, R., & Klauer, K. C. (2021). Partial derivatives for the
 * first-passage time distribution in Wiener diffusion models.
 * *Journal of Mathematical Psychology, 103*, 102550.
 * https://doi.org/10.1016/j.jmp.2021.102550
 * - Navarro, D. J., & Fuss, I. G. (2009). Fast and accurate calculations for
 * first-passage times in Wiener diffusion models.
 * *Journal of Mathematical Psychology, 53*(4), 222–230.
 * https://doi.org/10.1016/j.jmp.2009.02.003
 */

template <bool propto, typename T_y, typename T_alpha, typename T_delta,
          typename T_beta, typename T_t0, typename T_sv, typename T_sw,
          typename T_st0>
inline return_type_t<T_y, T_alpha, T_delta, T_beta, T_t0, T_sv, T_sw, T_st0>
wiener_full_lpdf(const T_y& y, const T_alpha& a, const T_delta& v,
                 const T_beta& w, const T_t0& t0, const T_sv& sv,
                 const T_sw& sw, const T_st0& st0) {
  double precision = 1e-4;

  return internal::wiener_full_prec_impl_lpdf<propto, T_y, T_alpha, T_delta,
                                              T_beta, T_t0, T_sv, T_sw, T_st0>(
      "wiener_full_lpdf", y, a, v, w, t0, sv, sw, st0, precision);
}

template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_t0, typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_alpha, T_delta, T_beta, T_t0, T_sv, T_sw, T_st0>
wiener_full_lpdf(const T_y& y, const T_alpha& a, const T_delta& v,
                 const T_beta& w, const T_t0& t0, const T_sv& sv,
                 const T_sw& sw, const T_st0& st0) {
  double precision = 1e-4;

  return internal::wiener_full_prec_impl_lpdf<false>(
      "wiener_full_lpdf", y, a, v, w, t0, sv, sw, st0, precision);
}

/** \ingroup prob_dists
 * The log of the first passage time density function for a (Wiener)
 * drift diffusion model with up to 7 parameters with the option
 * to control for the precision in the estimation of the partial derivatives.
 *
 * For \b Details see \c wiener_full_lpdf(). The usage and behavior of these
 * functions are the same excpet of the control over the precision.
 */

template <bool propto, typename T_y, typename T_alpha, typename T_delta,
          typename T_beta, typename T_t0, typename T_sv, typename T_sw,
          typename T_st0>
inline return_type_t<T_y, T_alpha, T_delta, T_beta, T_t0, T_sv, T_sw, T_st0>
wiener_full_prec_lpdf(const T_y& y, const T_alpha& a, const T_delta& v,
                      const T_beta& w, const T_t0& t0, const T_sv& sv,
                      const T_sw& sw, const T_st0& st0, const double& prec) {
  return internal::wiener_full_prec_impl_lpdf<propto, T_y, T_alpha, T_delta,
                                              T_beta, T_t0, T_sv, T_sw, T_st0>(
      "wiener_full_prec_lpdf", y, a, v, w, t0, sv, sw, st0, prec);
}

template <typename T_y, typename T_alpha, typename T_delta, typename T_beta,
          typename T_t0, typename T_sv, typename T_sw, typename T_st0>
inline return_type_t<T_y, T_alpha, T_delta, T_beta, T_t0, T_sv, T_sw, T_st0>
wiener_full_prec_lpdf(const T_y& y, const T_alpha& a, const T_delta& v,
                      const T_beta& w, const T_t0& t0, const T_sv& sv,
                      const T_sw& sw, const T_st0& st0, const double& prec) {
  return internal::wiener_full_prec_impl_lpdf<false>(
      "wiener_full_prec_lpdf", y, a, v, w, t0, sv, sw, st0, prec);
}

}  // namespace math
}  // namespace stan
#endif
