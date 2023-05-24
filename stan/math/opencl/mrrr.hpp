#ifndef STAN_MATH_OPENCL_MRRR_HPP
#define STAN_MATH_OPENCL_MRRR_HPP

#ifdef STAN_OPENCL

#include <cmath>
#include <queue>

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/prim/multiply.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/copy.hpp>
#include <stan/math/opencl/double_d.hpp>

#include <stan/math/opencl/kernels/mrrr.hpp>

namespace stan {
namespace math {
namespace internal {

using VectorXdd = Eigen::Matrix<double_d, -1, 1>;
using MatrixXdd = Eigen::Matrix<double_d, -1, -1>;

const double_d perturbation_range = 1e-20;

/**
 * Generates a random number for perturbing a relatively robust representation
 * @return A uniformly distributed random number between `1 - perturbation_range
 * / 2` and `1 + perturbation_range / 2`.
 */
inline double_d get_random_perturbation_multiplier() {
  static const double_d rand_norm = perturbation_range / RAND_MAX;
  static const double_d almost_one = 1 - perturbation_range * 0.5;
  return almost_one + std::rand() * rand_norm;
}

/**
 * Calculates bounds on eigenvalues of a symmetric tridiagonal matrix T using
 * Gresgorin discs.
 * @param diagonal Diagonal of T
 * @param subdiagonal Subdiagonal of T
 * @param[out] min_eigval Lower bound on eigenvalues.
 * @param[out] max_eigval Upper bound on eigenvalues.
 */
inline void get_gresgorin(const Eigen::Ref<const Eigen::VectorXd> diagonal,
                          const Eigen::Ref<const Eigen::VectorXd> subdiagonal,
                          double& min_eigval, double& max_eigval) {
  using std::fabs;
  const int n = diagonal.size();
  min_eigval = diagonal[0] - fabs(subdiagonal[0]);
  max_eigval = diagonal[0] + fabs(subdiagonal[0]);
  for (int i = 1; i < n - 1; i++) {
    min_eigval = std::min(min_eigval, diagonal[i] - fabs(subdiagonal[i])
                                          - fabs(subdiagonal[i - 1]));
    max_eigval = std::max(max_eigval, diagonal[i] + fabs(subdiagonal[i])
                                          + fabs(subdiagonal[i - 1]));
  }
  min_eigval = std::min(min_eigval, diagonal[n - 1] - fabs(subdiagonal[n - 2]));
  max_eigval = std::max(max_eigval, diagonal[n - 1] + fabs(subdiagonal[n - 2]));
}

/**
 * NaN-favouring max. If either of arguments is NaN, returns NaN. Otherwise
 * returns larger of the arguments.
 * @param a first argument
 * @param b second argument
 * @return NaN or larger argument
 */
inline double max_nan(double a, double b) { return isnan(a) || a > b ? a : b; }

/**
 * Calculates LDL decomposition of a shifted triagonal matrix T. D is diagonal,
 * L is lower unit triangular (diagonal elements are 1, all elements except
 * diagonal and subdiagonal are 0),T - shift * I = L * D * L^T. Also calculates
 * element growth of D: max(abs(D)).
 * @param diagonal Diagonal of T
 * @param subdiagonal Subdiagonal of T.
 * @param shift Shift.
 * @param[out] l Subdiagonal of L.
 * @param[out] d_plus Diagonal of D.
 * @return Element growth.
 */
inline double get_ldl(const Eigen::Ref<const Eigen::VectorXd> diagonal,
                      const Eigen::Ref<const Eigen::VectorXd> subdiagonal,
                      const double shift, VectorXdd& l, VectorXdd& d_plus) {
  using std::fabs;
  d_plus[0] = diagonal[0] - shift;
  double element_growth = fabs(d_plus[0].high);
  for (int i = 0; i < subdiagonal.size(); i++) {
    l[i] = subdiagonal[i] / d_plus[i];
    d_plus[i + 1] = diagonal[i + 1] - shift - l[i] * subdiagonal[i];
    element_growth = max_nan(element_growth, fabs(d_plus[i + 1].high));
  }
  return element_growth;
}

/**
 * Finds a good value for shift of the initial LDL factorization T - shift * I =
 * L * D * L^T.
 * @param diagonal Diagonal of T.
 * @param subdiagonal Subdiagonal of T.
 * @param l0 Subdiagonal of L.
 * @param d0 Diagonal of D.
 * @param min_eigval Lower bound on eigenvalues of T.
 * @param max_eigval High bound on eigenvalues of T
 * @param max_ele_growth Maximum desired element growth.
 * @return
 */
inline double find_initial_shift(
    const Eigen::Ref<const Eigen::VectorXd> diagonal,
    const Eigen::Ref<const Eigen::VectorXd> subdiagonal, VectorXdd& l0,
    VectorXdd& d0, const double min_eigval, const double max_eigval,
    const double max_ele_growth) {
  double shift = (max_eigval + min_eigval) * 0.5;
  double element_growth = get_ldl(diagonal, subdiagonal, shift, l0, d0);
  if (element_growth < max_ele_growth) {
    return shift;
  }
  double plus = (max_eigval - min_eigval) * 1e-15;
  while (!(element_growth < max_ele_growth)) {  // if condition was flipped it
                                                // would be wrong for the case
                                                // where element_growth is nan
    plus *= -2;
    element_growth = get_ldl(diagonal, subdiagonal, shift + plus, l0, d0);
  }
  return shift + plus;
}

/**
 * Calculates Sturm count of a LDL decomposition of a tridiagonal matrix -
 * number of eigenvalues larger or equal to shift. Uses stqds - calculation of
 * shifted LDL decomposition algorithm and counts number of positive elements in
 * D.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param shift Shift.
 * @return Sturm count.
 */
inline int get_sturm_count_ldl(const VectorXdd& l, const VectorXdd& d,
                               const double_d shift) {
  using std::isinf;
  const int n = l.size();
  double_d s = -shift;
  double_d d_plus;
  int count = 0;
  for (int i = 0; i < n; i++) {
    d_plus = s + d.coeff(i);
    count += d_plus >= 0;
    if (isinf(s)) {  // this happens if d_plus==0 -> in next
                     // iteration d_plus==inf and s==inf
      s = l.coeff(i) * l.coeff(i) * d.coeff(i) - shift;
    } else {
      s = l.coeff(i) * l.coeff(i) * s * (d.coeff(i) / d_plus) - shift;
    }
  }
  d_plus = s + d.coeff(n);
  count += d_plus >= 0;
  return count;
}

/**
 * Refines bounds on the i-th largest eigenvalue of LDL decomposition using
 * bisection.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param[in,out] low Low bound on the eigenvalue.
 * @param[in,out] high High bound on the eigenvalue.
 * @param i i-th eigenvalue
 */
inline void eigenval_bisect_refine(const VectorXdd& l, const VectorXdd& d,
                                   double_d& low, double_d& high, const int i) {
  using std::fabs;
  const double_d eps = 3e-20;
  while ((high - low) > eps * fabs(high + low)
         && fabs(high - low)
                > std::numeric_limits<double>::min()) {  // second term is for
                                                         // the case where the
                                                         // eigenvalue is 0 and
                                                         // division yields NaN
    double_d mid = (high + low) * 0.5;
    if (get_sturm_count_ldl(l, d, mid) > i) {
      low = mid;
    } else {
      high = mid;
    }
  }
}

/**
 * Shifts a LDL decomposition. The algorithm is sometimes called stationary
 * quotients-differences with shifts (stqds). D and D+ are diagonal, L and L+
 * are lower unit triangular (diagonal elements are 1, all elements except
 * diagonal and subdiagonal are 0). L * D * L^T - shift * I = L+ * D * L+^T.
 * Also calculates element growth of D+: max(abs(D+)).
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param shift Shift.
 * @param[out] l_plus Subdiagonal of L+.
 * @param[out] d_plus Diagonal of D+.
 * @return Element growth.
 */
inline double get_shifted_ldl(const VectorXdd& l, const VectorXdd& d,
                              const double_d shift, VectorXdd& l_plus,
                              VectorXdd& d_plus) {
  using std::fabs;
  using std::isinf;
  const int n = l.size();
  double_d s = -shift;
  double element_growth = 0;
  for (int i = 0; i < n; i++) {
    d_plus[i] = s + d[i];
    element_growth = max_nan(element_growth, fabs(d_plus[i].high));
    l_plus[i] = l[i] * (d[i] / d_plus[i]);
    if (isinf(d_plus[i]) && isinf(s)) {  // this happens if d_plus[i]==0 -> in
                                         // next iteration d_plus==inf and
                                         // s==inf
      s = l[i] * l[i] * d[i] - shift;
    } else {
      s = l_plus[i] * l[i] * s - shift;
    }
  }
  d_plus[n] = s + d[n];
  element_growth = max_nan(element_growth, fabs(d_plus[n].high));
  return element_growth;
}

/**
 * Finds good shift and shifts a LDL decomposition so as to keep element growth
 * low. L * D * L^T - shift * I = L2 * D2 * L2^T.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param low Low bound on wanted shift.
 * @param high High bound on wanted shift.
 * @param max_ele_growth Maximum desired element growth. If no better options
 * are found, it might be exceeded.
 * @param max_shift Maximal difference of shift from wanted bounds.
 * @param[out] l2 Subdiagonal of L2.
 * @param[out] d2 Diagonal of D2.
 * @param[out] shift Shift.
 * @param[out] min_element_growth Element growth achieved with resulting shift.
 */
inline void find_shift(const VectorXdd& l, const VectorXdd& d,
                       const double_d low, const double_d high,
                       const double max_ele_growth, const double_d max_shift,
                       VectorXdd& l2, VectorXdd& d2, double_d& shift,
                       double& min_element_growth) {
  VectorXdd l3(l2.size()), d3(d2.size());
  const std::vector<double_d> shifts = {
      low,
      high - max_shift * 0.1,
      low + max_shift * 0.1,
      high - max_shift * 0.25,
      low + max_shift * 0.25,
      high - max_shift * 0.5,
      low + max_shift * 0.5,
      high - max_shift * 0.75,
      low + max_shift * 0.75,
      high - max_shift,
      low + max_shift,
  };
  min_element_growth = std::numeric_limits<double>::infinity();
  for (double_d sh : shifts) {
    const double element_growth = get_shifted_ldl(l, d, sh, l3, d3);
    if (element_growth < min_element_growth) {
      l2.swap(l3);
      d2.swap(d3);
      shift = sh;
      min_element_growth = element_growth;
      if (element_growth <= max_ele_growth) {
        break;
      }
    }
  }
}

struct mrrr_task {
  int start, end;
  double_d shift;  // total shift, not just the last one
  VectorXdd l, d;
  int level;
};

/**
 * Calculates eigenvalues and eigenvectors of a irreducible tridiagonal matrix T
 * using multiple relatively robust representations (MRRR) algorithm. Use
 * `tridiagonal_eigensolver` if any subdiagonal element might be (very close to)
 * zero.
 * @param diagonal Diagonal of of T.
 * @param subdiagonal Subdiagonal of T.
 * @param[out] eigenvalues Eigenvlues.
 * @param[out] eigenvectors Eigenvectors.
 * @param min_rel_sep Minimal relative separation of eigenvalues before
 * computing eigenvectors.
 * @param maximum_ele_growth Maximal desired element growth of LDL
 * decompositions.
 */
template <bool need_eigenvectors = true>
inline void mrrr_cl(const Eigen::Ref<const Eigen::VectorXd> diagonal,
                    const Eigen::Ref<const Eigen::VectorXd> subdiagonal,
                    Eigen::Ref<Eigen::VectorXd> eigenvalues,
                    Eigen::Ref<Eigen::MatrixXd> eigenvectors,
                    const double min_rel_sep = 1e-4,
                    const double maximum_ele_growth = 15) {
  using std::copysign;
  using std::fabs;
  const double shift_error = 1e-19;
  const int n = diagonal.size();
  double min_eigval;
  double max_eigval;
  const Eigen::VectorXd subdiagonal_squared
      = subdiagonal.array() * subdiagonal.array();
  get_gresgorin(diagonal, subdiagonal, min_eigval, max_eigval);
  if (!need_eigenvectors) {
    matrix_cl<double> diagonal_cl(diagonal);
    matrix_cl<double> subdiagonal_squared_cl(subdiagonal_squared);
    matrix_cl<double> eigenvalues_cl(n, 1);
    matrix_cl<double_d> l_cl, d_cl, high_cl, low_cl(n, 1);
    opencl_kernels::eigenvals(cl::NDRange(n), diagonal_cl,
                              subdiagonal_squared_cl, l_cl, d_cl,
                              eigenvalues_cl, low_cl, high_cl, min_eigval,
                              max_eigval, 0, need_eigenvectors);
    eigenvalues = from_matrix_cl(eigenvalues_cl);
    return;
  }
  VectorXdd l(n - 1), d(n);
  const double max_ele_growth = maximum_ele_growth * (max_eigval - min_eigval);
  const double shift0 = find_initial_shift(
      diagonal, subdiagonal, l, d, min_eigval, max_eigval, max_ele_growth);
  for (int i = 0; i < n; i++) {
    if (i != n - 1) {
      l[i] = l[i] * get_random_perturbation_multiplier();
    }
    d[i] = d[i] * get_random_perturbation_multiplier();
  }
  VectorXdd high(n), low(n);

  matrix_cl<double> diagonal_cl(diagonal);
  matrix_cl<double> subdiagonal_squared_cl(subdiagonal_squared);
  matrix_cl<double_d> l_cl(l);
  matrix_cl<double_d> d_cl(d);
  matrix_cl<double> eigenvalues_cl(n, 1);
  matrix_cl<double_d> high_cl(n, 1);
  matrix_cl<double_d> low_cl(n, 1);
  opencl_kernels::eigenvals(cl::NDRange(n), diagonal_cl, subdiagonal_squared_cl,
                            l_cl, d_cl, eigenvalues_cl, low_cl, high_cl,
                            min_eigval, max_eigval, shift0, need_eigenvectors);
  eigenvalues = from_matrix_cl(eigenvalues_cl);
  high = from_matrix_cl(high_cl);
  low = from_matrix_cl(low_cl);

  MatrixXdd l_big(n - 1, n), d_big(n, n);

  std::queue<mrrr_task> block_queue;
  block_queue.push(mrrr_task{0, n, {shift0, 0}, std::move(l), std::move(d), 0});
  l.resize(n - 1);  // after move out
  d.resize(n);
  while (!block_queue.empty()) {
    const mrrr_task block = block_queue.front();
    block_queue.pop();
    double_d shift = std::numeric_limits<double>::infinity();
    double min_element_growth = std::numeric_limits<double>::infinity();
    VectorXdd l2(n - 1), d2(n), l_plus(n - 1), u_minus(n - 1);
    for (int i = block.start; i < block.end; i++) {
      int cluster_end;
      for (cluster_end = i + 1; cluster_end < block.end; cluster_end++) {
        const int prev = cluster_end - 1;
        const double_d end_threshold
            = low[prev] * (1 - copysign(min_rel_sep, low[prev]));
        if (high[cluster_end] < end_threshold) {
          break;
        }
      }
      cluster_end--;  // now this is the index of the last element of the
                      // cluster
      if (cluster_end > i) {  // cluster
        double_d a = high[cluster_end - 1], b = low[cluster_end];
        double_d max_shift
            = (high[cluster_end - 1] - low[cluster_end]) / min_rel_sep;
        double_d next_shift;
        double min_ele_growth;
        find_shift(block.l, block.d, low[cluster_end], high[i], max_ele_growth,
                   max_shift, l, d, next_shift, min_ele_growth);
        for (int j = i; j <= cluster_end; j++) {
          low[j] = low[j] * (double_d(1.0, 0.0) - copysign(shift_error, low[j]))
                   - next_shift - 1e-200;
          high[j]
              = high[j] * (double_d(1.0, 0.0) + copysign(shift_error, high[j]))
                - next_shift + 1e-200;
          eigenval_bisect_refine(l, d, low[j], high[j], j);
        }
        block_queue.push(mrrr_task{i, cluster_end + 1, block.shift + next_shift,
                                   std::move(l), std::move(d),
                                   block.level + 1});
        l.resize(n - 1);  // after move out
        d.resize(n);

        i = cluster_end;
      } else {  // isolated eigenvalue
        int twist_idx;
        const double_d low_gap
            = i == block.start
                  ? double_d(std::numeric_limits<double>::infinity())
                  : low[i - 1] - high[i];
        const double_d high_gap
            = i == block.end - 1
                  ? double_d(std::numeric_limits<double>::infinity())
                  : low[i] - high[i + 1];
        const double_d min_gap = std::min(low_gap, high_gap);
        const VectorXdd *l_ptr, *d_ptr;
        if (!(fabs(min_gap / ((high[i] + low[i]) * 0.5)) > min_rel_sep)) {
          if (!(fabs(min_gap / ((high[i] + low[i]) * 0.5 - shift)) > min_rel_sep
                && min_element_growth < max_ele_growth)) {
            const double_d max_shift = min_gap / min_rel_sep;
            find_shift(block.l, block.d, low[i], high[i], max_ele_growth,
                       max_shift, l2, d2, shift, min_element_growth);
          }
          low[i] = low[i] * (1 - copysign(shift_error, low[i])) - shift;
          high[i] = high[i] * (1 + copysign(shift_error, high[i])) - shift;
          eigenval_bisect_refine(l2, d2, low[i], high[i], i);
          l_ptr = &l2;
          d_ptr = &d2;
        } else {
          l_ptr = &block.l;
          d_ptr = &block.d;
        }
        l_big.col(i) = *l_ptr;
        d_big.col(i) = *d_ptr;
      }
    }
  }
  matrix_cl<double_d> l_big_cl(l_big.transpose());
  matrix_cl<double_d> d_big_cl(d_big.transpose());
  matrix_cl<double> subdiag_cl(subdiagonal);
  matrix_cl<double_d> shifted_eigvals_cl((low + high) * 0.5);
  matrix_cl<double_d> l_plus_cl(n - 1, n), u_minus_cl(n - 1, n), temp_cl(n, n);
  matrix_cl<double> eigenvectors_cl(n, n);
  opencl_kernels::get_eigenvectors(cl::NDRange(n), l_big_cl, d_big_cl,
                                   subdiag_cl, shifted_eigvals_cl, l_plus_cl,
                                   u_minus_cl, temp_cl, eigenvectors_cl);
  eigenvectors = from_matrix_cl(transpose((eigenvectors_cl)));
}

/**
 * Calculates eigenvalues and eigenvectors of a tridiagonal matrix T using MRRR
 * algorithm. If a subdiagonal element is close to zero compared to neighbors on
 * diagonal, the problem can be split into smaller ones.
 * @param diagonal_cl Diagonal of of T.
 * @param subdiagonal_cl Subdiagonal of T.
 * @param[out] eigenvalues_cl Eigenvlues.
 * @param[out] eigenvectors_cl Eigenvectors.
 * @param split_threshold Threshold for splitting the problem
 */
template <bool need_eigenvectors = true>
inline void tridiagonal_eigensolver_cl(const matrix_cl<double>& diagonal_cl,
                                       const matrix_cl<double>& subdiagonal_cl,
                                       matrix_cl<double>& eigenvalues_cl,
                                       matrix_cl<double>& eigenvectors_cl,
                                       const double split_threshold = 1e-15) {
  using std::fabs;
  using std::sqrt;
  const int n = diagonal_cl.rows();
  Eigen::MatrixXd eigenvectors;
  if (need_eigenvectors) {
    eigenvectors.resize(n, n);
  }
  Eigen::VectorXd eigenvalues(n);
  Eigen::VectorXd diagonal = from_matrix_cl<Eigen::VectorXd>(diagonal_cl);
  Eigen::VectorXd subdiagonal = from_matrix_cl<Eigen::VectorXd>(subdiagonal_cl);
  int last = 0;
  for (int i = 0; i < subdiagonal.size(); i++) {
    if (fabs(subdiagonal[i]) < split_threshold * sqrt(fabs(diagonal[i]))
                                   * sqrt(fabs(diagonal[i + 1]))) {
      if (need_eigenvectors) {
        eigenvectors.block(last, i + 1, i + 1 - last, n - i - 1)
            = Eigen::MatrixXd::Constant(i + 1 - last, n - i - 1, 0);
        eigenvectors.block(i + 1, last, n - i - 1, i + 1 - last)
            = Eigen::MatrixXd::Constant(n - i - 1, i + 1 - last, 0);
        if (last == i) {
          eigenvectors(last, last) = 1;
          eigenvalues[last] = diagonal[last];
        } else {
          mrrr_cl<need_eigenvectors>(
              diagonal.segment(last, i + 1 - last),
              subdiagonal.segment(last, i - last),
              eigenvalues.segment(last, i + 1 - last),
              eigenvectors.block(last, last, i + 1 - last, i + 1 - last));
        }
      } else {
        if (last == i) {
          eigenvalues[last] = diagonal[last];
        } else {
          mrrr_cl<need_eigenvectors>(diagonal.segment(last, i + 1 - last),
                                     subdiagonal.segment(last, i - last),
                                     eigenvalues.segment(last, i + 1 - last),
                                     eigenvectors);
        }
      }
      last = i + 1;
    }
  }
  if (last == n - 1) {
    if (need_eigenvectors) {
      eigenvectors(last, last) = 1;
    }
    eigenvalues[last] = diagonal[last];
  } else {
    if (need_eigenvectors) {
      mrrr_cl<need_eigenvectors>(
          diagonal.segment(last, n - last),
          subdiagonal.segment(last, subdiagonal.size() - last),
          eigenvalues.segment(last, n - last),
          eigenvectors.block(last, last, n - last, n - last));
    } else {
      mrrr_cl<need_eigenvectors>(
          diagonal.segment(last, n - last),
          subdiagonal.segment(last, subdiagonal.size() - last),
          eigenvalues.segment(last, n - last), eigenvectors);
    }
  }
  eigenvalues_cl = to_matrix_cl(std::move(eigenvalues));
  if (need_eigenvectors) {
    eigenvectors_cl = to_matrix_cl(std::move(eigenvectors));
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
#endif
