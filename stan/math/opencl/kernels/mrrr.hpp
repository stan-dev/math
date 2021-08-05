#ifndef STAN_MATH_GPU_KERNELS_MRRR_HPP
#define STAN_MATH_GPU_KERNELS_MRRR_HPP

#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/double_d.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
static const char* eigenvals_bisect_kernel_code = STRINGIFY(
    // \endcond

    /**
     * Calculates lower Sturm count of a tridiagonal matrix T - number of
     * eigenvalues lower than shift.
     * @param diagonal Diagonal of T.
     * @param subdiagonal_squared Squared elements of subdiagonal of T.
     * @param shift
     * @param n Size of T.
     * @return Sturm count
     */
    int get_sturm_count_tri(const __global double* diagonal,
                            const __global double* subdiagonal_squared,
                            const double shift, const int n) {
      double d = diagonal[0] - shift;
      int count = d >= 0.0;
      for (int j = 1; j < n; j++) {
        d = diagonal[j] - shift - subdiagonal_squared[j - 1] / d;
        count += d >= 0.0;
      }
      return count;
    }

    /**
     * Calculates i-th largest eigenvalue of tridiagonal matrix represented by a
     * LDL decomposition using bisection.
     * @param l Subdiagonal of L.
     * @param d Diagonal of D.
     * @param[out] low_res Resulting low bounds on eigenvalues.
     * @param[out] high_res Resulting high bounds on eigenvalues.
     * @param min_eigval Lower bound on all eigenvalues.
     * @param max_eigval Upper bound on all eigenvalues.
     * @param n Length of D.
     * @param i calculate i-th eigenvalue
     */
    void eigenvals_bisect(const __global double* diagonal,
                          const __global double* subdiagonal_squared,
                          double* low_res, double* high_res,
                          const double min_eigval, const double max_eigval,
                          const int n, const int i) {
      const double eps = 2 * DBL_EPSILON;

      double low = min_eigval;
      double high = max_eigval;

      while ((high - low) > eps * fabs(high + low)
             && fabs(high - low) > DBL_MIN) {
        double mid = (high + low) * 0.5;
        int count = get_sturm_count_tri(diagonal, subdiagonal_squared, mid, n);
        if (count > i) {
          low = mid;
        } else {
          high = mid;
        }
      }
      *low_res = low;
      *high_res = high;
    }

    /**
     * Calculates Sturm count of a LDL decomposition of a tridiagonal matrix -
     * number of eigenvalues larger or equal to shift. Uses stqds - calculation
     * of shifted LDL decomposition algorithm and counts number of positive
     * elements in D.
     * @param l Subdiagonal of L.
     * @param d Diagonal of D.
     * @param shift Shift.
     * @param n Length of L.
     * @return Sturm count.
     */
    int get_sturm_count_ldl(const __global double_d* l,
                            const __global double_d* d, const double_d shift,
                            const int n) {
      double_d s = neg(shift);
      double_d l_plus;
      double_d d_plus;
      int count = 0;
      for (int i = 0; i < n; i++) {
        d_plus = add_dd_dd(s, d[i]);
        count += ge_dd_d(d_plus, 0.0);
        if (isinf_dd(s)) {  // this happens if d_plus==0 -> in next iteration
                            // d_plus==inf and s==inf
          s = sub_dd_dd(mul_dd_dd(mul_dd_dd(l[i], l[i]), d[i]), shift);
        } else {
          s = sub_dd_dd(mul_dd_dd(mul_dd_dd(mul_dd_dd(l[i], l[i]), s),
                                  div_dd_dd(d[i], d_plus)),
                        shift);
        }
      }
      d_plus = add_dd_dd(s, d[n]);
      count += ge_dd_d(d_plus, 0);
      return count;
    }

    /**
     * Refines bounds on the i-th largest eigenvalue of a LDL decomposition
     * using bisection.
     * @param l Subdiagonal of L.
     * @param d Diagonal of D.
     * @param[in,out] low_res Low bound on the eigenvalue.
     * @param[in,out] high_res High bound on the eigenvalue.
     * @param n size of `d`
     * @param i i-th eigenvalue
     */
    void eigenvals_bisect_refine(const __global double_d* l,
                                 const __global double_d* d, double_d* low_res,
                                 double_d* high_res, const int n, const int i) {
      double_d eps = (double_d){3e-20, 0};
      double_d min_norm = (double_d){DBL_MIN, 0};

      double_d low = *low_res;
      double_d high = *high_res;

      while (gt_dd_dd(
                 abs_dd(div_dd_dd(sub_dd_dd(high, low), add_dd_dd(high, low))),
                 eps)
             && gt_dd_dd(abs_dd(sub_dd_dd(high, low)), min_norm)) {
        double_d mid = mul_dd_d(add_dd_dd(high, low), 0.5);
        int count = get_sturm_count_ldl(l, d, mid, n - 1);
        if (count > i) {
          low = mid;
        } else {
          high = mid;
        }
      }
      *low_res = low;
      *high_res = high;
    }

    /**
     * Calculates eigenvalues of a tridiagonal matrix T and refines shifted
     * eigenvalues using shifted LDL decomposition of T.
     * @param diagonal diagonal of T
     * @param subdiagonal_squared element-wise squared subdiagonal of T
     * @param l subdiagonal of L
     * @param d diagonal of D
     * @param[out] eigval_global eigenvalues of T
     * @param[out] shifted_low_global lower bounds on shifted eigenvalues
     * @param[out] shifted_high_global upper bounds on shifted eigenvalues
     * @param min_eigval initial lower bound on eigenvalues
     * @param max_eigval initial upper bound on eigenvalues
     * @param shift shift of the LDL decomposition
     */
    __kernel void eigenvals(
        const __global double* diagonal,
        const __global double* subdiagonal_squared, const __global double_d* l,
        const __global double_d* d, __global double* eigval_global,
        __global double_d* shifted_low_global,
        __global double_d* shifted_high_global, const double min_eigval,
        const double max_eigval, const double shift, const char do_refine) {
      const int i = get_global_id(0);
      const int n = get_global_size(0);

      double low_eig, high_eig;
      eigenvals_bisect(diagonal, subdiagonal_squared, &low_eig, &high_eig,
                       min_eigval, max_eigval, n, i);
      eigval_global[i] = (low_eig + high_eig) * 0.5;
      if (do_refine) {
        double_d low_shifted = (double_d){low_eig - shift, 0};
        double_d high_shifted = (double_d){high_eig - shift, 0};
        low_shifted = mul_dd_dd(
            low_shifted,
            sub_dd_d((double_d){1, 0}, copysign_d_dd(1e-18 * n, low_shifted)));
        high_shifted = mul_dd_dd(
            high_shifted,
            add_dd_d((double_d){1, 0}, copysign_d_dd(1e-18 * n, high_shifted)));
        eigenvals_bisect_refine(l, d, &low_shifted, &high_shifted, n, i);
        shifted_low_global[i] = low_shifted;
        shifted_high_global[i] = high_shifted;
      }
    }
    // \cond
);
// \endcond

const kernel_cl<in_buffer, in_buffer, in_buffer, in_buffer, out_buffer,
                out_buffer, out_buffer, double, double, double, char>
    eigenvals("eigenvals", {stan::math::internal::double_d_src,
                            eigenvals_bisect_kernel_code});

// \cond
static const char* get_eigenvectors_kernel_code = STRINGIFY(
    // \endcond

    /**
     * Calculates shifted LDL and UDU factorizations. Combined with twist index
     * they form twisted factorization for calculation of an eigenvector
     * corresponding to eigenvalue that is equal to the shift. Tha algorithm is
     * sometimes called diferential twisted quotient-differences with shifts
     * (dtwqds). L * D * L^T - shift * I = L+ * D+ * L+^T = U- * D- * U-^T D, D+
     * and D- are diagonal, L and L+ are lower unit triangular (diagonal
     * elements are 1, all elements except diagonal and subdiagonal are 0), U-
     * is upper unit triangular (diagonal elements are 1, all elements except
     * diagonal and superdiagonal are 0)
     * @param l Subdiagonal of L.
     * @param d Diagonal of D.
     * @param shift Shift.
     * @param[out] l_plus Subdiagonal of L+.
     * @param[out] u_minus Superdiagonal of U-.
     * @param s Temporary array of the same size as d.
     * @return Twist index.
     */
    int get_twisted_factorization(
        const __global double_d* l, const __global double_d* d, double_d shift,
        __global double_d* l_plus, __global double_d* u_minus,
        __global double_d* s) {
      int n = get_global_size(0);
      int gid = get_global_id(0);
      int m = n - 1;
      // calculate shifted ldl
      s[gid] = neg(shift);
      for (int i = 0; i < m; i++) {
        double_d d_plus = add_dd_dd(s[i * n + gid], d[i * n + gid]);
        l_plus[i * n + gid]
            = mul_dd_dd(l[i * n + gid], div_dd_dd(d[i * n + gid], d_plus));
        if (isnan_dd(l_plus[i * n + gid])) {  // d_plus==0
          // one (or both) of d[i], l[i] is very close to 0
          if (lt_dd_dd(abs_dd(l[i * n + gid]), abs_dd(d[i * n + gid]))) {
            l_plus[i * n + gid]
                = mul_dd_d(d[i * n + gid], copysign_d_dd(1., l[i * n + gid])
                                               * copysign_d_dd(1., d_plus));
          } else {
            l_plus[i * n + gid]
                = mul_dd_d(l[i * n + gid], copysign_d_dd(1., d[i * n + gid])
                                               * copysign_d_dd(1., d_plus));
          }
        }
        s[(i + 1) * n + gid] = sub_dd_dd(
            mul_dd_dd(mul_dd_dd(l_plus[i * n + gid], l[i * n + gid]),
                      s[i * n + gid]),
            shift);
        if (isnan_dd(s[(i + 1) * n + gid])) {
          if (gt_dd_dd(abs_dd(l_plus[i * n + gid]),
                       abs_dd(s[i * n + gid]))) {  // l_plus[i * n + gid] == inf
            if (gt_dd_dd(abs_dd(s[i * n + gid]),
                         abs_dd(l[i * n + gid]))) {  // l[i*n+gid]==0
              s[(i + 1) * n + gid] = sub_dd_dd(
                  mul_dd_d(s[i * n + gid],
                           copysign_d_dd(1., l[i * n + gid])
                               * copysign_d_dd(1., l_plus[i * n + gid])),
                  shift);
            } else {  // s[i*n+gid]==0
              s[(i + 1) * n + gid] = sub_dd_dd(
                  mul_dd_d(l[i * n + gid],
                           copysign_d_dd(1., s[i * n + gid])
                               * copysign_d_dd(1., l_plus[i * n + gid])),
                  shift);
            }
          } else {  // s[i*n+gid]==inf
            if (gt_dd_dd(abs_dd(l_plus[i * n + gid]),
                         abs_dd(l[i * n + gid]))) {  // l[i]==0
              s[(i + 1) * n + gid]
                  = sub_dd_dd(mul_dd_d(l_plus[i * n + gid],
                                       copysign_d_dd(1., l[i * n + gid])
                                           * copysign_d_dd(1., s[i * n + gid])),
                              shift);
            } else {  // l_plus[i]==0
              s[(i + 1) * n + gid] = sub_dd_dd(
                  mul_dd_d(l[i * n + gid],
                           copysign_d_dd(1., s[i * n + gid])
                               * copysign_d_dd(1., l_plus[i * n + gid])),
                  shift);
            }
          }
        }
      }
      // calculate shifted udu and twist index
      double_d p = sub_dd_dd(d[m * n + gid], shift);
      double_d min_gamma = abs_dd(add_dd_dd(s[m * n + gid], d[m * n + gid]));
      int twist_index = m;

      for (int i = m - 1; i >= 0; i--) {
        double_d d_minus
            = add_dd_dd(mul_dd_dd(mul_dd_dd(d[i * n + gid], l[i * n + gid]),
                                  l[i * n + gid]),
                        p);
        double_d t = div_dd_dd(d[i * n + gid], d_minus);
        u_minus[i * n + gid] = mul_dd_dd(l[i * n + gid], t);
        if (isnan_dd(u_minus[i * n + gid])) {
          if (isnan_dd(t)) {
            double t_high = copysign_d_dd(1., d[i * n + gid])
                            * copysign_d_dd(1., d_minus);
            t.high = t_high;
            t.low = 0;
            u_minus[i * n + gid] = mul_dd_d(l[i * n + gid], t_high);
          } else {  // t==inf, l[i*n+gid]==0
            u_minus[i * n + gid]
                = mul_dd_d(d[i * n + gid], copysign_d_dd(1., l[i * n + gid])
                                               * copysign_d_dd(1., t));
          }
        }
        double_d gamma = abs_dd(add_dd_dd(s[i * n + gid], mul_dd_dd(t, p)));
        if (isnan_dd(gamma)) {  // t==inf, p==0 OR t==0, p==inf
          double_d d_sign
              = mul_dd_d(d[i * n + gid],
                         copysign_d_dd(1., d_minus) * copysign_d_dd(1., t));
          gamma = abs_dd(add_dd_dd(s[i * n + gid], d_sign));
          p = sub_dd_dd(d_sign, shift);
        } else {  // usual case
          p = sub_dd_dd(mul_dd_dd(p, t), shift);
        }
        if (lt_dd_dd(gamma, min_gamma)) {
          min_gamma = gamma;
          twist_index = i;
        }
      }
      return twist_index;
    }

    /**
     * Calculates an eigenvector from twisted factorization T - shift * I = L+
     * * D+ * L+^T = U- * D- * U-^T.
     * @param l_plus Subdiagonal of the L+.
     * @param u_minus Superdiagonal of the U-.
     * @param subdiag Subdiagonal of T
     * @param twist_idx Twist index.
     * @param[out] eigenvectors Matrix in which to store resulting vectors.
     */
    void calculate_eigenvector(const __global double_d* l_plus,
                               const __global double_d* u_minus,
                               const __global double* subdiag, int twist_idx,
                               __global double* eigenvectors) {
      int n = get_global_size(0);
      int gid = get_global_id(0);
      int i = gid;
      eigenvectors[twist_idx * n + gid] = 1;
      double_d last = (double_d){1, 0};
      double_d last2 = (double_d){1, 0};
      double norm = 1;
      // part of the eigenvector after the twist index
      for (int j = twist_idx + 1; j < n; j++) {
        if (last.high != 0 || last.low != 0) {
          last2 = last;
          last = neg(mul_dd_dd(u_minus[(j - 1) * n + gid], last));
          eigenvectors[j * n + gid] = last.high;
        } else {
          double_d tmp = last;
          last = mul_dd_d(last2, -subdiag[j - 2] / subdiag[j - 1]);
          last2 = tmp;
          if (isnan(last.high) || isinf(last.high)) {  // subdiag[j - 1]==0
            last = (double_d){0, 0};
          }
          eigenvectors[j * n + gid] = last.high;
        }
        norm += eigenvectors[j * n + gid] * eigenvectors[j * n + gid];
      }
      last = (double_d){eigenvectors[twist_idx * n + gid], 0};
      last2 = (double_d){1, 0};
      // part of the eigenvector before the twist index
      for (int j = twist_idx - 1; j >= 0; j--) {
        if (last.high != 0 || last.low != 0) {
          last2 = last;
          last = neg(mul_dd_dd(l_plus[j * n + gid], last));
          eigenvectors[j * n + gid] = last.high;
        } else {
          double_d tmp = last;
          last = mul_dd_d(last2, -subdiag[j + 1] / subdiag[j]);
          if (isnan(last.high) || isinf(last.high)) {  // subdiag[j]==0
            last = (double_d){0, 0};
          }
          eigenvectors[j * n + gid] = last.high;
        }
        norm += eigenvectors[j * n + gid] * eigenvectors[j * n + gid];
      }
      norm = 1 / sqrt(norm);
      // normalize the eigenvector
      for (int j = 0; j < n; j++) {
        eigenvectors[j * n + gid] *= norm;
      }
    }

    /**
     * Calculates eigenvectors for (shifted) eigenvalues.
     * @param l Each row is a subdiagonal of the matrix L from LDL
     * decomposition for one eigenvalue.
     * @param d Each row is a diagonal of the matrix D from LDL
     * decomposition for one eigenvalue.
     * @param subdiag Subdiagonal of the tridiagonal matrix.
     * @param shifted_eigvals shifted eigenvalues
     * @param l_plus Temporary array of the same size as l
     * @param u_minus Temporary array of the same size as l
     * @param temp Temporary array of the same size as d
     * @param[out] eigenvectors Each row is one eigenvector.
     */
    __kernel void get_eigenvectors(
        const __global double_d* l, const __global double_d* d,
        const __global double* subdiag,
        const __global double_d* shifted_eigvals, __global double_d* l_plus,
        __global double_d* u_minus, __global double_d* temp,
        __global double* eigenvectors) {
      int twist_idx = get_twisted_factorization(
          l, d, shifted_eigvals[get_global_id(0)], l_plus, u_minus, temp);
      calculate_eigenvector(l_plus, u_minus, subdiag, twist_idx, eigenvectors);
    }
    // \cond
);
// \endcond

const kernel_cl<in_buffer, in_buffer, in_buffer, in_buffer, in_out_buffer,
                in_out_buffer, in_out_buffer, out_buffer>
    get_eigenvectors("get_eigenvectors", {stan::math::internal::double_d_src,
                                          get_eigenvectors_kernel_code});

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif
