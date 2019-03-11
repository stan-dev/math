#ifndef STAN_MATH_GPU_KERNELS_MRRR_HPP
#define STAN_MATH_GPU_KERNELS_MRRR_HPP
#ifndef STAN_OPENCL
#error "NO STAN_OPENCL"
#endif
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_cl.hpp>

namespace stan {
namespace math {
namespace opencl_kernels {

// \cond
const char* eigenvals_bisect_kernel_code = STRINGIFY(
// \endcond
        /**
         * Calculates Sturm count of a LDL decomposition of a tridiagonal matrix - number of eigenvalues larger or equal to shift.
         * Uses stqds - calculation of shifted LDL decomposition algorithm and counts number of positive elements in D.
         * @param l Subdiagonal of L.
         * @param d Diagonal of D.
         * @param shift Shift.
         * @return Sturm count.
         */
        int getSturmCountLdl(__global double* l, const __global double* d, double shift, int n){
          double s = -shift;
          double l_plus;
          double d_plus;
          int count = 0;
          for (int i = 0; i < n; i++) {
            d_plus = s + d[i];
            count += d_plus >= 0;
            if (isinf(d_plus) && isinf(s)) { // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
              s = l[i] * l[i] * d[i] - shift;
            }
            else {
              s = l[i] * l[i] * s * (d[i] / d_plus) - shift;
            }
          }
          d_plus = s + d[n];
          count += d_plus >= 0;
          return count;
        }

        /**
         * Calculates eigenvalues of tridiagonal matrix represented by a LDL decomposition using bisection.
         * @param l Subdiagonal of L.
         * @param d Diagonal of D.
         * @param[out] low_res Resulting low bounds on eigenvalues.
         * @param[out] high_res Resulting high bounds on eigenvalues.
         * @param min_eigval Lower bound on all eigenvalues.
         * @param max_eigval Upper bound on all eigenvalues.
         */
        __kernel void eigenvals_bisect(__global double* l, const __global double* d, __global double* low_res, __global double* high_res,
                const double min_eigval, const double max_eigval, const int n) {
          const int lid = get_local_id(0);
          const int gid = get_global_id(0);
          const int gsize = get_global_size(0);
          const int lsize = get_local_size(0);
          const int ngroups = get_num_groups(0);
          const int wgid = get_group_id(0);

          double eps = 3e-16;
          double min_norm = 3e-308; //(approximately) smallest normalized double, larger than 0

          int i=gid;
          double low = min_eigval;
          double high = max_eigval;

          while (fabs((high - low) / (high + low)) > eps && fabs(high - low) > min_norm) {
            double mid = (high + low) * 0.5;
            int count = getSturmCountLdl(l, d, mid, n-1);
            if (count > i) {
              low = mid;
            }
            else {
              high = mid;
            }
          }
          low_res[i]=low;
          high_res[i]=high;
        }
// \cond
);
// \endcond

// \cond
const char* get_eigenvectors_kernel_code = STRINGIFY(
// \endcond

        /**
         * Shifts a LDL decomposition. The algorithm is sometimes called stationary quotients-differences with shifts (stqds).
         * D and D+ are diagonal, L and L+ are lower unit triangular (diagonal elements are 1,
         * all elements except diagonal and subdiagonal are 0). L * D * L^T - shift * I = L+ * D * L+^T.
         * Also calculates element growth of D+: sum(abs(D+)) / abs(sum(D+)).
         * @param l Subdiagonal of L.
         * @param d Diagonal of D.
         * @param shift Shift.
         * @param[out] l_plus Subdiagonal of L+.
         * @param[out] d_plus Diagonal of D+.
         * @return Element growth.
         */
        double get_shifted_ldl(const __global double* l, const __global double* d, double shift, __global double* l_plus, __global double* d_plus) {
          int n=get_global_size(0);
          int m=n-1;
          int gid=get_global_id(0);
          double s = -shift;
          double element_growth = 0;
          double element_growth_denominator = 0;
          for (int i = 0; i < m; i++) {
            d_plus[i*n+gid] = s + d[i*n+gid];
            element_growth += fabs(d_plus[i*n+gid]);
            element_growth_denominator += d_plus[i*n+gid];
            l_plus[i*n+gid] = l[i*n+gid] * (d[i*n+gid] / d_plus[i*n+gid]);
            if (isinf(d_plus[i*n+gid]) && isinf(s)) { // this happens if d_plus[i]==0 -> in next iteration d_plus==inf and s==inf
              s = l[i*n+gid] * l[i*n+gid] * d[i*n+gid] - shift;
            }
            else {
              s = l_plus[i*n+gid] * l[i*n+gid] * s - shift;
            }
          }
          d_plus[m*n+gid] = s + d[m*n+gid];
          element_growth += fabs(d_plus[m*n+gid]);
          return element_growth / fabs(element_growth_denominator);
        }

        /**
         * Swaps two pointers.
         * @param a First pointer.
         * @param b Second pointer.
         */
        void swap(__global double** a, __global double** b){
          __global double* c=*a;
          *a=*b;
          *b=c;
        }


        /**
         * Finds good shift and shifts a LDL decomposition so as to keep element growth low. L * D * L^T - shift * I = L2 * D2 * L2^T.
         * @param l Subdiagonal of L.
         * @param d Diagonal of D.
         * @param low Low bound on wanted shift.
         * @param high High bound on wanted shift.
         * @param max_ele_growth Maximum desired element growth. If no better options are found it might be exceeded.
         * @param max_shift Maximal difference of shhift from wanted bounds.
         * @param[out] l2 Subdiagonal of L2.
         * @param[out] d2 Diagonal of D2.
         * @param l3 Temporary array of the same size as d.
         * @param d3 Temporary array of the same size as d
         * @param[out] shift Shift.
         * @param[out] min_element_growth Element growth achieved with resulting shift.
         */
        void findShift(const __global double* l, const __global double* d, double low, double high, double max_ele_growth, double max_shift,
                       __global double** l2, __global double** d2, __global double** l3, __global double** d3, double* shift, double* min_element_growth
                       ) {
          double shifts[11];
          shifts[0] = low;
          shifts[1] = high - max_shift * 0.1;
          shifts[2] = low + max_shift * 0.1;
          shifts[3] = high - max_shift * 0.25;
          shifts[4] = low + max_shift * 0.25;
          shifts[5] = high - max_shift * 0.5;
          shifts[6] = low + max_shift * 0.5;
          shifts[7] = high - max_shift * 0.75;
          shifts[8] = low + max_shift * 0.75;
          shifts[9] = high - max_shift;
          shifts[10] = low + max_shift;

          *min_element_growth = INFINITY;
          for (int i=0;i<sizeof(shifts)/(sizeof(shifts[0]));i++) {
            //sh -= shift0;
            double element_growth = get_shifted_ldl(l, d, shifts[i], *l3, *d3);
            //cout << " element growth: " << element_growth << " at " << sh << endl;
            if (element_growth < *min_element_growth) {
              swap(l2,l3);
              swap(d2,d3);
              *shift = shifts[i];
              *min_element_growth = element_growth;
              if (element_growth <= max_ele_growth) {
                break;
              }
            }
          }
//  cout << "\t\t" << " element growth: " << min_element_growth << " at " << shift << endl;
        }

        /**
         * Calculates Sturm count of a LDL decomposition of a tridiagonal matrix - number of eigenvalues larger or equal to shift.
         * Uses stqds - calculation of shifted LDL decomposition algorithm and counts number of positive elements in D.
         * @param l Subdiagonal of L.
         * @param d Diagonal of D.
         * @param shift Shift.
         * @return Sturm count.
         */
        int getSturmCountLdl(const __global double* l, const __global double* d, double shift){
          int gid = get_global_id(0);
          int n = get_global_size(0);
          int m = n - 1;
          double s = -shift;
          double d_plus;
          int count = 0;
          for (int i = 0; i < m; i++) {
            d_plus = s + d[i*n+gid];
            count += d_plus >= 0;
            if (isinf(d_plus) && isinf(s)) { // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
              s = l[i*n+gid] * l[i*n+gid] * d[i*n+gid] - shift;
            }
            else {
              s = l[i*n+gid] * l[i*n+gid] * s * (d[i*n+gid] / d_plus) - shift;
            }
          }
          d_plus = s + d[m*n+gid];
          count += d_plus >= 0;
          return count;
        }

        /**
         * Refines bounds on eigenvalues of LDL decomposition of a matrix using bisection.
         * @param l Subdiagonal of L.
         * @param d Diagonal of D.
         * @param low[in,out] Low bound on the eigenvalue.
         * @param high[in,out] High bound on the eigenvalue.
         */
        void eigenvalBisectRefine(const __global double* l, const __global double* d, double* low, double* high) {
          int i=get_global_id(0);
          double eps = 3e-16;
          while (fabs((*high - *low) / (*high + *low)) > eps && fabs(*high - *low) > DBL_MIN) { // second term is for the case where the eigenvalue is 0 and division yields NaN
            double mid = (*high + *low) * 0.5;
            if (getSturmCountLdl(l, d, mid) > i) {
              *low = mid;
            }
            else {
              *high = mid;
            }
          }
        }

        /**
         * Calculates shifted LDL and UDU factorizations. Combined with twist index they form twisted factorization for calculation
         * of an eigenvector corresponding to eigenvalue that is equal to the shift. Tha algorithm is sometimes called diferential twisted quotient-differences with shifts (dtwqds).
         * L * D * L^T - shift * I = L+ * D+ * L+^T = U- * D- * U-^T
         * D, D+ and D- are diagonal, L and L+ are lower unit triangular (diagonal elements are 1, all elements except diagonal and subdiagonal are 0),
         * U- is upper unit triangular (diagonal elements are 1, all elements except diagonal and superdiagonal are 0)
         * @param l Subdiagonal of L.
         * @param d Diagonal of D.
         * @param shift Shift.
         * @param[out] l_plus Subdiagonal of L+.
         * @param[out] u_minus Superdiagonal of U-.
         * @param s Temporary array of the same size as d.
         * @return Twist index.
         */
        int get_twisted_factorization(const __global double* l, const __global double* d, double shift, __global double* l_plus, __global double* u_minus, __global double* s) {
          int n = get_global_size(0);
          int gid = get_global_id(0);
          int m = n - 1;
          //calculate shifted ldl
          s[gid] = -shift;
          for (int i = 0; i < m; i++) {
            double d_plus = s[i*n+gid] + d[i*n+gid];
            l_plus[i*n+gid] = l[i*n+gid] * (d[i*n+gid] / d_plus);
            if (isnan(l_plus[i*n+gid])) { //d_plus==0
              //one (or both) of d[i], l[i] is very close to 0
              if (fabs(l[i*n+gid]) < fabs(d[i*n+gid])) {
                l_plus[i*n+gid] = d[i*n+gid] * copysign(1., l[i*n+gid]) * copysign(1., d_plus);
              }
              else {
                l_plus[i*n+gid] = l[i*n+gid] * copysign(1., d[i*n+gid]) * copysign(1., d_plus);
              }
            }
//            printf("%d set l_plus to %lf", gid, l_plus[i*n+gid]);
            s[(i + 1)*n+gid] = l_plus[i*n+gid] * l[i*n+gid] * s[i*n+gid] - shift;
            if (isnan(s[(i + 1)*n+gid])) {
              if (fabs(l_plus[i*n+gid]) > fabs(s[i*n+gid])) { //l_plus[i*n+gid]==inf
                if (fabs(s[i*n+gid]) > fabs(l[i*n+gid])) { //l[i*n+gid]==0
                  s[(i + 1)*n+gid] = s[i*n+gid] * copysign(1., l[i*n+gid]) * copysign(1., l_plus[i*n+gid]) - shift;
                }
                else { //s[i*n+gid]==0
                  s[(i + 1)*n+gid] = l[i*n+gid] * copysign(1., s[i*n+gid]) * copysign(1., l_plus[i*n+gid]) - shift;
                }
              }
              else { //s[i*n+gid]==inf
                if (fabs(l_plus[i*n+gid]) > fabs(l[i*n+gid])) { //l[i]==0
                  s[(i + 1)*n+gid] = l_plus[i*n+gid] * copysign(1., l[i*n+gid]) * copysign(1., s[i*n+gid]) - shift;
                }
                else { //l_plus[i]==0
                  s[(i + 1)*n+gid] = l[i*n+gid] * copysign(1., s[i*n+gid]) * copysign(1., l_plus[i*n+gid]) - shift;
                }
              }
            }
//            printf("%d set s to %lf", gid, s[(i + 1)*n+gid]);
          }
          //calculate shifted udu and twist index
          double p = d[m*n+gid] - shift;
          double min_gamma = fabs(s[m*n+gid] + d[m*n+gid]);
          int twist_index = m;

          for (int i = m - 1; i >= 0; i--) {
            double d_minus = d[i*n+gid] * l[i*n+gid] * l[i*n+gid] + p;
            double t = d[i*n+gid] / d_minus;
            u_minus[i*n+gid] = l[i*n+gid] * t;
            if (isnan(u_minus[i*n+gid])) {
              if (isnan(t)) {
                t = copysign(1., d[i*n+gid]) * copysign(1., d_minus);
                u_minus[i*n+gid] = l[i*n+gid] * t;
              }
              else { //t==inf, l[i*n+gid]==0
                u_minus[i*n+gid] = d[i*n+gid] * copysign(1., l[i*n+gid]) * copysign(1., t);
              }
            }
//            printf("%d set u_minus to %lf", gid, u_minus[i*n+gid]);
            double gamma = fabs(s[i*n+gid] + t * p);
            if (isnan(gamma)) { //t==inf, p==0 OR t==0, p==inf
              double d_sign = d[i*n+gid] * copysign(1., d_minus) * copysign(1., t);
              gamma = fabs(s[i*n+gid] + d_sign);
              p = d_sign - shift;
            }
            else { //usual case
              p = p * t - shift;
            }
            if (gamma < min_gamma) {
              min_gamma = gamma;
              twist_index = i;
            }
          }
          return twist_index;
        }


        /**
          * Calculates eigenvectors from twisted factorization T - shift * I = L+ * D+ * L+^T = U- * D- * U-^T.
          * @param l_plus Subdiagonal of the L+.
          * @param u_minus Superdiagonal of the U-.
          * @param subdiag Subdiagonal of T
          * @param twist_idx Twist index.
          * @param[out] eigenvecs Matrix in which to store resulting vectors.
          */
        void calculateEigenvector(const __global double* l_plus, const __global double* u_minus, const __global double* subdiag, int twist_idx, __global double* eigenvecs) {
          int n = get_global_size(0);
          int gid = get_global_id(0);
          int i = gid;
          eigenvecs[twist_idx*n+gid] = 1;
          double norm=1;
          for (int j = twist_idx + 1; j < n; j++) {
            if (eigenvecs[(j - 1)*n+gid] != 0) {
              eigenvecs[j*n+gid] = -u_minus[(j - 1)*n+gid] * eigenvecs[(j - 1)*n+gid];
            }
            else {
              eigenvecs[j*n+gid] = -subdiag[j - 2] * eigenvecs[(j - 2)*n+gid] / subdiag[j - 1];
              if (isnan(eigenvecs[j*n+gid]) || isinf(eigenvecs[j*n+gid])) { //subdiag[j - 1]==0
                eigenvecs[j*n+gid] = 0;
              }
            }
            norm += eigenvecs[j*n+gid] * eigenvecs[j*n+gid];
          }
          for (int j = twist_idx - 1; j >= 0; j--) {
            if (eigenvecs[(j + 1)*n+gid] != 0) {
              eigenvecs[j*n+gid] = -l_plus[j*n+gid] * eigenvecs[(j + 1)*n+gid];
            }
            else {
              eigenvecs[j*n+gid] = -subdiag[j + 1] * eigenvecs[(j + 2)*n+gid] / subdiag[j];
              if (isnan(eigenvecs[j*n+gid]) || isinf(eigenvecs[j*n+gid])) { //subdiag[j]==0
                eigenvecs[j*n+gid] = 0;
              }
            }
            norm += eigenvecs[j*n+gid] * eigenvecs[j*n+gid];
          }
//          norm=1/sqrt(norm);
//          for(int j=0;j<n;j++){
//            eigenvecs[j*n+gid]*=norm;
//          }
        }

        /**
         * Calculates eigenvectors for distinct (shifted) eigenvalues.
         * @param subdiag Subdiagonal of the tridiagonal matrix.
         * @param l Each row is a subdiagonal of the matrix L from LDL decomposition for one eigenvalue.
         * @param d Each row is a diagonal of the matrix D from LDL decomposition for one eigenvalue.
         * @param low_glob Lower bounds on shifted eigenvalues
         * @param high_glob High bounds on shifted eigenvalues.
         * @param min_gap_glob Minimal absolute gap of each eigenvalue to the closest eigenvalue
         * @param l2 Temporary array of the same size as d
         * @param d2 Temporary array of the same size as d
         * @param temp1 Temporary array of the same size as d
         * @param temp2 Temporary array of the same size as d
         * @param temp3 Temporary array of the same size as d
         * @param eigenvecs Each row is one eigenvector.
         * @param min_rel_sep Minimal relative separation of eigenvalues before computing eigenvectors.
         * @param max_ele_growth Maximal desired element growth of LDL decompositions.
         */
        __kernel void get_eigenvectors(const __global double* subdiag, const __global double* l, const __global double* d,
                const __global double* low_glob, const __global double* high_glob, const __global double* min_gap_glob,
                __global double* l2, __global double* d2, __global double* temp1, __global double* temp2, __global double* temp3, __global double* eigenvecs,
                double min_rel_sep, double max_ele_growth) {
          const int gid = get_global_id(0);
          const int n = get_global_size(0);

          double shift_error = 1e-14;
          double min_gap = min_gap_glob[gid];
          double high = high_glob[gid];
          double low = low_glob[gid];


          const __global double* l_ptr;
          const __global double* d_ptr;
          if (!(fabs(min_gap / ((high + low) * 0.5)) > min_rel_sep)) {
            double max_shift = min_gap / min_rel_sep;
            double shift;
            double min_element_growth;
            findShift(l, d, low, high, max_ele_growth, max_shift, &l2, &d2, &temp1, &temp2, &shift, &min_element_growth);
//            printf("%d shift: %lf\n", gid, shift);
            low = low * (1 - copysign(shift_error, low)) - shift;
            high = high * (1 + copysign(shift_error, high)) - shift;
            eigenvalBisectRefine(l2, d2, &low, &high);
            l_ptr = l2;
            d_ptr = d2;
          }
          else {
//            printf("%d default ok\n", gid);
            l_ptr = l;
            d_ptr = d;
          }
//          printf("%d eigval: %lf %lf\n", gid, low, high);
          __global double* l_plus = temp1;
          __global double* u_minus = temp2;
//          for(int i=0;i<4;i++){
//            barrier(CLK_LOCAL_MEM_FENCE);
//            if(gid==i){
//              printf("%d: ",gid);
//              for(int j=0;j<3;j++){
//                printf("%lf, ",l_ptr[j*n+gid]);
//              }
//              printf("\n%d: ", gid);
//              for(int j=0;j<4;j++){
//                printf("%lf, ",d_ptr[j*n+gid]);
//              }
//              printf("\n");
//            }
//          }
          int twist_idx = get_twisted_factorization(l_ptr, d_ptr, (low + high) * 0.5, l_plus, u_minus, temp3);
//          for(int i=0;i<4;i++){
//            barrier(CLK_LOCAL_MEM_FENCE);
//            if(gid==i){
//              printf("%d: ",gid);
//              for(int j=0;j<3;j++){
//                printf("%lf, ",l_plus[j*n+gid]);
//              }
//              printf("\n%d: ", gid);
//              for(int j=0;j<3;j++){
//                printf("%lf, ",u_minus[j*n+gid]);
//              }
//              printf("\n");
//            }
//          }
          calculateEigenvector(l_plus, u_minus, subdiag, twist_idx, eigenvecs);
        }
// \cond
);
// \endcond

const global_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, double, double, int >
        eigenvals_bisect("eigenvals_bisect", eigenvals_bisect_kernel_code);

const global_range_kernel<cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, double, double>
        get_eigenvectors("get_eigenvectors", get_eigenvectors_kernel_code);

}  // namespace opencl_kernels
}  // namespace math
}  // namespace stan
#endif
#endif