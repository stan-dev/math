#ifndef STAN_MATH_GPU_KERNELS_MRRR_HPP
#define STAN_MATH_GPU_KERNELS_MRRR_HPP

#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>
#include <stan/math/opencl/double_d.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
const char* eigenvals_bisect_kernel_code = STRINGIFY(
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
     * Calculates eigenvalues of tridiagonal matrix represented by a LDL
     * decomposition using bisection.
     * @param l Subdiagonal of L.
     * @param d Diagonal of D.
     * @param[out] low_res Resulting low bounds on eigenvalues.
     * @param[out] high_res Resulting high bounds on eigenvalues.
     * @param min_eigval Lower bound on all eigenvalues.
     * @param max_eigval Upper bound on all eigenvalues.
     * @param n Length of D.
     * @param i calculate i-th eigenvalue
     */
    void eigenvals_bisect(
        __global double* diagonal, const __global double* subdiagonal_squared,
        double* low_res, double* high_res, const double min_eigval,
        const double max_eigval, const int n, const int i) {
      const double eps = 3e-16;
      const double min_norm
          = 3e-308;  //(approximately) smallest normalized double, larger than 0

      double low = min_eigval;
      double high = max_eigval;

      while (fabs((high - low) / (high + low)) > eps
             && fabs(high - low) > min_norm) {
        double mid = (high + low) * 0.5;
        int count = get_sturm_count_tri(diagonal, subdiagonal_squared, mid, n);
        //        printf("%d shift %lf count %d\n", i, mid, count);
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
     * @return Sturm count.
     * @param n Length of L.
     */
    int get_sturm_count_ldl(__global double_d* l, const __global double_d* d,
                            double_d shift, int n) {
      double_d s = neg(shift);
      double_d l_plus;
      double_d d_plus;
      int count = 0;
      for (int i = 0; i < n; i++) {
        d_plus = add_dd_dd(s, d[i]);
        count += ge_dd_d(d_plus, 0);
        if (isinf_dd(d_plus)
            && isinf_dd(s)) {  // this happens if d_plus==0 -> in next iteration
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

    void eigenvals_bisect_refine(__global double_d* l,
                                 const __global double_d* d, double_d* low_res,
                                 double_d* high_res, const int n, const int i) {
      double_d eps = (double_d){3e-16, 0};
      double_d min_norm = (double_d){
          3e-308,
          0};  //(approximately) smallest normalized double, larger than 0

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

    __kernel void eigenvals(
        __global double* diagonal, const __global double* subdiagonal_squared,
        __global double_d* l, const __global double_d* d,
        __global double* eigval_low_global, __global double* eigval_high_global,
        __global double_d* shifted_low_global,
        __global double_d* shifted_high_global, const double min_eigval,
        const double max_eigval, const double shift) {
      const int i = get_global_id(0);
      const int n = get_global_size(0);

      double low_eig, high_eig;
      eigenvals_bisect(diagonal, subdiagonal_squared, &low_eig, &high_eig,
                       min_eigval, max_eigval, n, i);
      //      printf("%d: %lf %lf\n", i, low_eig, high_eig);
      eigval_low_global[i] = low_eig;
      eigval_high_global[i] = high_eig;
      double_d low_shifted = (double_d){low_eig - shift, 0};
      double_d high_shifted = (double_d){high_eig - shift, 0};
      low_shifted
          = mul_dd_d(low_shifted, (1 - copysign_d_dd(1e-14 * n, low_shifted)));
      high_shifted = mul_dd_d(high_shifted,
                              (1 + copysign_d_dd(1e-14 * n, high_shifted)));
      eigenvals_bisect_refine(l, d, &low_shifted, &high_shifted, n, i);
      shifted_low_global[i] = low_shifted;
      shifted_high_global[i] = high_shifted;
    }
    // \cond
);
// \endcond

const kernel_cl<in_buffer, in_buffer, in_buffer, in_buffer, out_buffer,
                out_buffer, out_buffer, out_buffer, double, double, double>
    eigenvals("eigenvals", {stan::math::internal::double_d_src,
                            eigenvals_bisect_kernel_code});

// \cond
const char* get_eigenvectors_kernel_code = STRINGIFY(
    // \endcond
    /**
     * Shifts a LDL decomposition. The algorithm is sometimes called stationary
     * quotients-differences with shifts (stqds). D and D+ are diagonal, L and
     * L+ are lower unit triangular (diagonal elements are 1, all elements
     * except diagonal and subdiagonal are 0). L * D * L^T - shift * I = L+ * D
     * * L+^T. Also calculates element growth of D+: sum(abs(D+)) /
     * abs(sum(D+)).
     * @param l Subdiagonal of L.
     * @param d Diagonal of D.
     * @param shift Shift.
     * @param[out] l_plus Subdiagonal of L+.
     * @param[out] d_plus Diagonal of D+.
     * @return Element growth.
     */
    double get_shifted_ldl(const __global double* l, const __global double* d,
                           double shift, __global double* l_plus,
                           __global double* d_plus) {
      int n = get_global_size(0);
      int m = n - 1;
      int gid = get_global_id(0);
      double s = -shift;
      double element_growth = 0;
      double element_growth_denominator = 0;
      for (int i = 0; i < m; i++) {
        d_plus[i * n + gid] = s + d[i * n + gid];
        element_growth += fabs(d_plus[i * n + gid]);
        element_growth_denominator += d_plus[i * n + gid];
        l_plus[i * n + gid]
            = l[i * n + gid] * (d[i * n + gid] / d_plus[i * n + gid]);
        if (isinf(d_plus[i * n + gid])
            && isinf(s)) {  // this happens if d_plus[i]==0 -> in next iteration
                            // d_plus==inf and s==inf
          s = l[i * n + gid] * l[i * n + gid] * d[i * n + gid] - shift;
        } else {
          s = l_plus[i * n + gid] * l[i * n + gid] * s - shift;
        }
      }
      d_plus[m * n + gid] = s + d[m * n + gid];
      element_growth += fabs(d_plus[m * n + gid]);
      return element_growth / fabs(element_growth_denominator);
    }

    /**
     * Swaps two pointers.
     * @param a First pointer.
     * @param b Second pointer.
     */
    void swap(__global double** a, __global double** b) {
      __global double* c = *a;
      *a = *b;
      *b = c;
    }

    //    /**
    //     * Finds good shift and shifts a LDL decomposition so as to keep
    //     element
    //     * growth low. L * D * L^T - shift * I = L2 * D2 * L2^T.
    //     * @param l Subdiagonal of L.
    //     * @param d Diagonal of D.
    //     * @param low Low bound on wanted shift.
    //     * @param high High bound on wanted shift.
    //     * @param max_ele_growth Maximum desired element growth. If no better
    //     * options are found it might be exceeded.
    //     * @param max_shift Maximal difference of shhift from wanted bounds.
    //     * @param[out] l2 Subdiagonal of L2.
    //     * @param[out] d2 Diagonal of D2.
    //     * @param l3 Temporary array of the same size as d.
    //     * @param d3 Temporary array of the same size as d
    //     * @param[out] shift Shift.
    //     * @param[out] min_element_growth Element growth achieved with
    //     resulting
    //     * shift.
    //     */
    //    void find_shift(
    //        const __global double* l, const __global double* d, double low,
    //        double high, double max_ele_growth, double max_shift,
    //        __global double** l2, __global double** d2, __global double** l3,
    //        __global double** d3, double* shift, double* min_element_growth) {
    //      double shifts[11];
    //      shifts[0] = low;
    //      shifts[1] = high - max_shift * 0.1;
    //      shifts[2] = low + max_shift * 0.1;
    //      shifts[3] = high - max_shift * 0.25;
    //      shifts[4] = low + max_shift * 0.25;
    //      shifts[5] = high - max_shift * 0.5;
    //      shifts[6] = low + max_shift * 0.5;
    //      shifts[7] = high - max_shift * 0.75;
    //      shifts[8] = low + max_shift * 0.75;
    //      shifts[9] = high - max_shift;
    //      shifts[10] = low + max_shift;

    //      *min_element_growth = INFINITY;
    //      for (int i = 0; i < sizeof(shifts) / (sizeof(shifts[0])); i++) {
    //        double element_growth = get_shifted_ldl(l, d, shifts[i], *l3,
    //        *d3); if (element_growth < *min_element_growth) {
    //          swap(l2, l3);
    //          swap(d2, d3);
    //          *shift = shifts[i];
    //          *min_element_growth = element_growth;
    //          if (element_growth <= max_ele_growth) {
    //            break;
    //          }
    //        }
    //      }
    //    }

    //    /**
    //     * Calculates Sturm count of a LDL decomposition of a tridiagonal
    //     matrix -
    //     * number of eigenvalues larger or equal to shift. Uses stqds -
    //     calculation
    //     * of shifted LDL decomposition algorithm and counts number of
    //     positive
    //     * elements in D.
    //     * @param l Subdiagonal of L.
    //     * @param d Diagonal of D.
    //     * @param shift Shift.
    //     * @return Sturm count.
    //     */
    //    int get_sturm_count_ldl(const __global double* l, const __global
    //    double* d,
    //                            double shift) {
    //      int gid = get_global_id(0);
    //      int n = get_global_size(0);
    //      int m = n - 1;
    //      double s = -shift;
    //      double d_plus;
    //      int count = 0;
    //      for (int i = 0; i < m; i++) {
    //        d_plus = s + d[i * n + gid];
    //        count += d_plus >= 0;
    //        if (isinf(d_plus)
    //            && isinf(s)) {  // this happens if d_plus==0 -> in next
    //                            // iteration d_plus==inf and s==inf
    //          s = l[i * n + gid] * l[i * n + gid] * d[i * n + gid] - shift;
    //        } else {
    //          s = l[i * n + gid] * l[i * n + gid] * s * (d[i * n + gid] /
    //          d_plus)
    //              - shift;
    //        }
    //      }
    //      d_plus = s + d[m * n + gid];
    //      count += d_plus >= 0;
    //      return count;
    //    }

    //    /**
    //     * Refines bounds on eigenvalues of LDL decomposition of a matrix
    //     using
    //     * bisection.
    //     * @param l Subdiagonal of L.
    //     * @param d Diagonal of D.
    //     * @param[in,out] low Low bound on the eigenvalue.
    //     * @param[in,out] high High bound on the eigenvalue.
    //     */
    //    void eigenval_bisect_refine(const __global double* l,
    //                                const __global double* d, double* low,
    //                                double* high) {
    //      int i = get_global_id(0);
    //      double eps = 3e-16;
    //      while (fabs((*high - *low) / (*high + *low)) > eps
    //             && fabs(*high - *low)
    //                    > DBL_MIN) {  // second term is for the case where the
    //                                  // eigenvalue is 0 and division yields
    //                                  NaN
    //        double mid = (*high + *low) * 0.5;
    //        if (get_sturm_count_ldl(l, d, mid) > i) {
    //          *low = mid;
    //        } else {
    //          *high = mid;
    //        }
    //      }
    //    }

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
//        if (isnan_dd(l_plus[i * n + gid])) {  // d_plus==0
//          // one (or both) of d[i], l[i] is very close to 0
//          if (lt_dd_dd(abs_dd(l[i * n + gid]), abs_dd(d[i * n + gid]))) {
//            l_plus[i * n + gid]
//                = mul_dd_d(d[i * n + gid], copysign_d_dd(1., l[i * n + gid])
//                                               * copysign_d_dd(1., d_plus));
//          } else {
//            l_plus[i * n + gid]
//                = mul_dd_d(l[i * n + gid], copysign_d_dd(1., d[i * n + gid])
//                                               * copysign_d_dd(1., d_plus));
//          }
//        }
        s[(i + 1) * n + gid] = sub_dd_dd(
            mul_dd_dd(mul_dd_dd(l_plus[i * n + gid], l[i * n + gid]),
                      s[i * n + gid]),
            shift);
//        if (isnan_dd(s[(i + 1) * n + gid])) {
//          if (gt_dd_dd(abs_dd(l_plus[i * n + gid]),
//                       abs_dd(s[i * n + gid]))) {  // l_plus[i*n+gid]==inf
//            if (gt_dd_dd(abs_dd(s[i * n + gid]),
//                         abs_dd(l[i * n + gid]))) {  // l[i*n+gid]==0
//              s[(i + 1) * n + gid] = sub_dd_dd(
//                  mul_dd_d(s[i * n + gid],
//                           copysign_d_dd(1., l[i * n + gid])
//                               * copysign_d_dd(1., l_plus[i * n + gid])),
//                  shift);
//            } else {  // s[i*n+gid]==0
//              s[(i + 1) * n + gid] = sub_dd_dd(
//                  mul_dd_d(l[i * n + gid],
//                           copysign_d_dd(1., s[i * n + gid])
//                               * copysign_d_dd(1., l_plus[i * n + gid])),
//                  shift);
//            }
//          } else {  // s[i*n+gid]==inf
//            if (gt_dd_dd(abs_dd(l_plus[i * n + gid]),
//                         abs_dd(l[i * n + gid]))) {  // l[i]==0
//              s[(i + 1) * n + gid]
//                  = sub_dd_dd(mul_dd_d(l_plus[i * n + gid],
//                                       copysign_d_dd(1., l[i * n + gid])
//                                           * copysign_d_dd(1., s[i * n + gid])),
//                              shift);
//            } else {  // l_plus[i]==0
//              s[(i + 1) * n + gid] = sub_dd_dd(
//                  mul_dd_d(l[i * n + gid],
//                           copysign_d_dd(1., s[i * n + gid])
//                               * copysign_d_dd(1., l_plus[i * n + gid])),
//                  shift);
//            }
//          }
//        }
      }
      // calculate shifted udu and twist index
      double_d p = sub_dd_dd(d[m * n + gid], shift);
      double_d min_gamma = abs_dd(add_dd_dd(s[m * n + gid], d[m * n + gid]));
//      printf("GPU %d gamma %le\n", gid, min_gamma.high);
      int twist_index = m;

      for (int i = m - 1; i >= 0; i--) {
        double_d d_minus
            = add_dd_dd(mul_dd_dd(mul_dd_dd(d[i * n + gid], l[i * n + gid]),
                                  l[i * n + gid]),
                        p);
        double_d t = div_dd_dd(d[i * n + gid], d_minus);
        u_minus[i * n + gid] = mul_dd_dd(l[i * n + gid], t);
//        if (isnan_dd(u_minus[i * n + gid])) {
//          if (isnan_dd(t)) {
//            double t_high = copysign_d_dd(1., d[i * n + gid])
//                            * copysign_d_dd(1., d_minus);
//            t.high = t_high;
//            t.low = 0;
//            u_minus[i * n + gid] = mul_dd_d(l[i * n + gid], t_high);
//          } else {  // t==inf, l[i*n+gid]==0
//            u_minus[i * n + gid]
//                = mul_dd_d(d[i * n + gid], copysign_d_dd(1., l[i * n + gid])
//                                               * copysign_d_dd(1., t));
//          }
//        }
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
//                    printf("GPU %d %le\n", gid, s[(i + 1) * n + gid].high);
//        printf("GPU %d gamma %le\n", gid, gamma.high);
        if (lt_dd_dd(gamma, min_gamma)) {
          min_gamma = gamma;
          twist_index = i;
        }
      }
      //                if(gid==0) printf("GPU %d %d\n", twist_index);
      return twist_index;
    }

    /**
     * Calculates eigenvectors from twisted factorization T - shift * I = L+
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
      for (int j = 0; j < n; j++) {
        eigenvectors[j * n + gid] *= norm;
      }
    }

    /**
     * Calculates eigenvectors for distinct (shifted) eigenvalues.
     * @param l Each row is a subdiagonal of the matrix L from LDL
     * decomposition for one eigenvalue.
     * @param d Each row is a diagonal of the matrix D from LDL
     * decomposition for one eigenvalue.
     * @param subdiag Subdiagonal of the tridiagonal matrix.
     * @param shifted_eigvals shifted eigenvalues
     * @param temp Temporary array of the same size as d
     * @param eigenvectors Each row is one eigenvector.
     * decompositions.
     */
    __kernel void get_eigenvectors(
        const __global double_d* l, const __global double_d* d,
        const __global double* subdiag,
        const __global double_d* shifted_eigvals, __global double_d* l_plus,
        __global double_d* u_minus, __global double_d* temp,
        __global double* eigenvectors) {
      int twist_idx = get_twisted_factorization(
          l, d, shifted_eigvals[get_global_id(0)], l_plus, u_minus, temp);
//      printf("GPU %d %d\n", (int)get_global_id(0), twist_idx);
      //      for (int i = 0; i < get_global_size(0) - 1; i++) {
      //        printf("%d %le %le\n", (int)get_global_id(0),
      //               l_plus[i * get_global_size(0) + get_global_id(0)].high,
      //               u_minus[i * get_global_size(0) + get_global_id(0)].high);
      //      }
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
