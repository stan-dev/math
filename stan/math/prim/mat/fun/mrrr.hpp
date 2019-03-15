#ifndef STAN_MATH_PRIM_MAT_FUN_MRRR_HPP
#define STAN_MATH_PRIM_MAT_FUN_MRRR_HPP

#include <queue>
#include <cmath>

#include <Eigen/Dense>

namespace stan {
namespace math {
namespace internal {

const double perturbation_range = 1e-15;

/**
 * Generates a random number for perturbing a relatively robust representation
 * @return A uniformly distributed random number between `1 - perturbation_range
 * / 2` and `1 + perturbation_range / 2`.
 */
inline double get_random_perturbation_multiplier() {
  static const double rand_norm = perturbation_range / RAND_MAX;
  static const double almost_one = 1 - perturbation_range * 0.5;
  return almost_one + std::rand() * rand_norm;
}

/**
 * Calculates LDL decomposition of a shifted triagonal matrix T. D is diagonal,
 * L is lower unit triangular (diagonal elements are 1, all elements except
 * diagonal and subdiagonal are 0),T - shift * I = L * D * L^T. Also calculates
 * element growth of D: sum(abs(D)) / abs(sum(D)).
 * @param diagonal Diagonal of T
 * @param subdiagonal Subdiagonal of T.
 * @param shift Shift.
 * @param[out] l Subdiagonal of L.
 * @param[out] d_plus Diagonal of D.
 * @return Element growth.
 */
double get_ldl(const Eigen::Ref<const Eigen::VectorXd> diagonal,
               const Eigen::Ref<const Eigen::VectorXd> subdiagonal,
               const double shift, Eigen::VectorXd& l,
               Eigen::VectorXd& d_plus) {
  using std::fabs;
  d_plus[0] = diagonal[0] - shift;
  double element_growth = fabs(d_plus[0]);
  double element_growth_denominator = d_plus[0];
  for (int i = 0; i < subdiagonal.size(); i++) {
    l[i] = subdiagonal[i] / d_plus[i];
    d_plus[i] *= get_random_perturbation_multiplier();
    d_plus[i + 1] = diagonal[i + 1] - shift - l[i] * subdiagonal[i];
    l[i] *= get_random_perturbation_multiplier();
    element_growth += fabs(d_plus[i + 1]);
    element_growth_denominator += d_plus[i + 1];
  }
  d_plus[subdiagonal.size()] *= get_random_perturbation_multiplier();
  return element_growth / fabs(element_growth_denominator);
}

/**
 * Shifts a LDL decomposition. The algorithm is sometimes called stationary
 * quotients-differences with shifts (stqds). D and D+ are diagonal, L and L+
 * are lower unit triangular (diagonal elements are 1, all elements except
 * diagonal and subdiagonal are 0). L * D * L^T - shift * I = L+ * D * L+^T.
 * Also calculates element growth of D+: sum(abs(D+)) / abs(sum(D+)).
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param shift Shift.
 * @param[out] l_plus Subdiagonal of L+.
 * @param[out] d_plus Diagonal of D+.
 * @return Element growth.
 */
double get_shifted_ldl(const Eigen::VectorXd& l, const Eigen::VectorXd& d,
                       const double shift, Eigen::VectorXd& l_plus,
                       Eigen::VectorXd& d_plus) {
  using std::fabs;
  using std::isinf;
  const int n = l.size();
  double s = -shift;
  double element_growth = 0;
  double element_growth_denominator = 0;
  for (int i = 0; i < n; i++) {
    d_plus[i] = s + d[i];
    element_growth += fabs(d_plus[i]);
    element_growth_denominator += d_plus[i];
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
  element_growth += fabs(d_plus[n]);
  return element_growth / fabs(element_growth_denominator);
}

/**
 * Calculates shifted LDL and UDU factorizations. Combined with twist index they
 * form twisted factorization for calculation of an eigenvector corresponding to
 * eigenvalue that is equal to the shift. Tha algorithm is sometimes called
 * diferential twisted quotient-differences with shifts (dtwqds). L * D * L^T -
 * shift * I = L+ * D+ * L+^T = U- * D- * U-^T D, D+ and D- are diagonal, L and
 * L+ are lower unit triangular (diagonal elements are 1, all elements except
 * diagonal and subdiagonal are 0), U- is upper unit triangular (diagonal
 * elements are 1, all elements except diagonal and superdiagonal are 0)
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param shift Shift.
 * @param[out] l_plus Subdiagonal of L+.
 * @param[out] u_minus Superdiagonal of U-.
 * @return Twist index.
 */
int get_twisted_factorization(const Eigen::VectorXd& l,
                              const Eigen::VectorXd& d, const double shift,
                              Eigen::VectorXd& l_plus,
                              Eigen::VectorXd& u_minus) {
  using std::copysign;
  using std::fabs;
  using std::isnan;
  const int n = l.size();
  // calculate shifted ldl
  Eigen::VectorXd s(n + 1);
  s[0] = -shift;
  for (int i = 0; i < n; i++) {
    double d_plus = s[i] + d[i];
    l_plus[i] = l[i] * (d[i] / d_plus);
    if (isnan(l_plus[i])) {  // d_plus==0
      if (fabs(l[i])
          < fabs(d[i])) {  // one (or both) of d[i], l[i] is very close to 0
        l_plus[i] = d[i] * copysign(1., l[i]) * copysign(1., d_plus);
      } else {
        l_plus[i] = l[i] * copysign(1., d[i]) * copysign(1., d_plus);
      }
    }
    s[i + 1] = l_plus[i] * l[i] * s[i] - shift;
    if (isnan(s[i + 1])) {
      if (fabs(l_plus[i]) > fabs(s[i])) {  // l_plus[i]==inf
        if (fabs(s[i]) > fabs(l[i])) {     // l[i]==0
          s[i + 1]
              = s[i] * copysign(1., l[i]) * copysign(1., l_plus[i]) - shift;
        } else {  // s[i]==0
          s[i + 1]
              = l[i] * copysign(1., s[i]) * copysign(1., l_plus[i]) - shift;
        }
      } else {                               // s[i]==inf
        if (fabs(l_plus[i]) > fabs(l[i])) {  // l[i]==0
          s[i + 1]
              = l_plus[i] * copysign(1., l[i]) * copysign(1., s[i]) - shift;
        } else {  // l_plus[i]==0
          s[i + 1]
              = l[i] * copysign(1., s[i]) * copysign(1., l_plus[i]) - shift;
        }
      }
    }
  }
  // calculate shifted udu and twist index
  double p = d[n] - shift;
  double min_gamma = fabs(s[n] + d[n]);
  int twist_index = n;

  for (int i = n - 1; i >= 0; i--) {
    double d_minus = d[i] * l[i] * l[i] + p;
    double t = d[i] / d_minus;
    u_minus[i] = l[i] * t;
    if (isnan(u_minus[i])) {
      if (isnan(t)) {
        t = copysign(1., d[i]) * copysign(1., d_minus);
        u_minus[i] = l[i] * t;
      } else {  // t==inf, l[i]==0
        u_minus[i] = d[i] * copysign(1., l[i]) * copysign(1., t);
      }
    }
    double gamma = fabs(s[i] + t * p);
    if (isnan(gamma)) {  // t==inf, p==0 OR t==0, p==inf
      double d_sign = d[i] * copysign(1., d_minus) * copysign(1., t);
      gamma = fabs(s[i] + d_sign);
      p = d_sign - shift;
    } else {  // general case
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
 * Calculates Sturm count of a LDL decomposition of a tridiagonal matrix -
 * number of eigenvalues larger or equal to shift. Uses stqds - calculation of
 * shifted LDL decomposition algorithm and counts number of positive elements in
 * D.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param shift Shift.
 * @return Sturm count.
 */
int get_sturm_count_ldl(const Eigen::VectorXd& l, const Eigen::VectorXd& d,
                        const double shift) {
  using std::isinf;
  const int n = l.size();
  double s = -shift;
  double d_plus;
  int count = 0;
  for (int i = 0; i < n; i++) {
    d_plus = s + d[i];
    count += d_plus >= 0;
    if (isinf(d_plus) && isinf(s)) {  // this happens if d_plus==0 -> in next
                                      // iteration d_plus==inf and s==inf
      s = l[i] * l[i] * d[i] - shift;
    } else {
      s = l[i] * l[i] * s * (d[i] / d_plus) - shift;
    }
  }
  d_plus = s + d[n];
  count += d_plus >= 0;
  return count;
}

/**
 * Refines bounds on the i-th largest eigenvalue of LDL decomposition of a
 * matrix using bisection.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param low[in,out] Low bound on the eigenvalue.
 * @param high[in,out] High bound on the eigenvalue.
 * @param i i-th eigenvalue
 */
void eigenval_bisect_refine(const Eigen::VectorXd& l, const Eigen::VectorXd& d,
                            double& low, double& high, const int i) {
  using std::fabs;
  const double eps = 3e-16;
  while (fabs((high - low) / (high + low)) > eps
         && fabs(high - low)
                > std::numeric_limits<double>::min()) {  // second term is for
                                                         // the case where the
                                                         // eigenvalue is 0 and
                                                         // division yields NaN
    double mid = (high + low) * 0.5;
    if (get_sturm_count_ldl(l, d, mid) > i) {
      low = mid;
    } else {
      high = mid;
    }
  }
}

/**
 * Calculates bounds on eigenvalues of a symmetric tridiagonal matrix T using
 * Gresgorin discs.
 * @param diagonal Diagonal of T
 * @param subdiagonal Subdiagonal of T
 * @param min_eigval[out] Lower bound on eigenvalues.
 * @param max_eigval[out] Upper bound on eigenvalues.
 */
void get_gresgorin(const Eigen::Ref<const Eigen::VectorXd> diagonal,
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

const int BISECT_K = 8;

/**
 * Calculates lower Sturm count of a tridiagonal matrix T - number of
 * eigenvalues lower than shift for up to BISECT_K different shifts.
 * @param diagonal Diagonal of T.
 * @param subdiagonal_squared Squared elements of subdiagonal of T.
 * @param shifts Up to 8 different shifts. First `n_valid` are used.
 * @param n_valid How many Sturm counts to actually compute.
 * @return Array of Sturm counts of size BISECT_K. First `n_valid` are actual
 * results.
 */
Eigen::Array<int, BISECT_K, 1> get_sturm_count_T_vec(
    const Eigen::Ref<const Eigen::VectorXd> diagonal,
    const Eigen::VectorXd& subdiagonal_squared,
    const Eigen::Array<double, BISECT_K, 1>& shifts, const int n_valid) {
  Eigen::Array<double, BISECT_K, 1> d;
  d.head(n_valid) = diagonal[0] - shifts.head(n_valid);
  Eigen::Array<int, BISECT_K, 1> counts;
  counts.head(n_valid) = (d.head(n_valid) < 0).cast<int>();
  for (int j = 1; j < diagonal.size(); j++) {
    d.head(n_valid) = diagonal[j] - shifts.head(n_valid)
                      - subdiagonal_squared[j - 1] / d.head(n_valid);
    counts.head(n_valid) += (d.head(n_valid) < 0).cast<int>();
  }
  return counts;
}

struct bisection_task {
  int start, end;
  double low, high;
};

/**
 * Calculates eigenvalues of tridiagonal matrix T using bisection.
 * @param diagonal Diagonal of T.
 * @param subdiagonal_squared Squared elements of the subdiagonal.
 * @param min_eigval Lower bound on all eigenvalues.
 * @param max_eigval Upper bound on all eigenvalues.
 * @param low[out] Lower bounds on eigenvalues.
 * @param high[out] Upper bounds on eigenvalues.
 */
void eigenvals_bisect(const Eigen::Ref<const Eigen::VectorXd> diagonal,
                      const Eigen::VectorXd& subdiagonal_squared,
                      const double min_eigval, const double max_eigval,
                      Eigen::VectorXd& low, Eigen::VectorXd& high) {
  using std::fabs;
  const int n = diagonal.size();
  const double eps = 3e-16;

  std::queue<bisection_task> task_queue;
  task_queue.push(bisection_task{0, n, min_eigval, max_eigval});
  while (!task_queue.empty()) {
    const int n_valid = std::min(BISECT_K, static_cast<int>(task_queue.size()));
    Eigen::Array<double, BISECT_K, 1> shifts;
    bisection_task t[BISECT_K];
    for (int i = 0; i < n_valid; i++) {
      t[i] = task_queue.front();
      task_queue.pop();
    }
    for (int i = 0; i < BISECT_K; i++) {
      const int task_idx = i % n_valid;
      const int idx_in_task = i / n_valid;
      const int task_total
          = BISECT_K / n_valid + (BISECT_K % n_valid > task_idx);
      shifts[i] = t[task_idx].low
                  + (t[task_idx].high - t[task_idx].low) * (idx_in_task + 1.)
                        / (task_total + 1);
    }
    const Eigen::Array<int, BISECT_K, 1> counts = get_sturm_count_T_vec(
        diagonal, subdiagonal_squared, shifts, BISECT_K);
    for (int i = 0; i < n_valid; i++) {
      if (counts[i] >= t[i].start + 1) {
        if ((t[i].high - shifts[i]) / fabs(shifts[i]) > eps
            && shifts[i] - t[i].low > std::numeric_limits<double>::min()) {
          task_queue.push({t[i].start, counts[i], t[i].low, shifts[i]});
        } else {
          const int n_eq = counts[i] - t[i].start;
          low.segment(t[i].start, n_eq)
              = Eigen::VectorXd::Constant(n_eq, t[i].low);
          high.segment(t[i].start, n_eq)
              = Eigen::VectorXd::Constant(n_eq, shifts[i]);
        }
      }
    }
    for (int i = 0; i < BISECT_K; i++) {
      const int task_idx = i % n_valid;
      const int idx_in_task = i / n_valid;
      const int task_total
          = BISECT_K / n_valid + (BISECT_K % n_valid > task_idx);
      int my_end = t[task_idx].end;
      double my_high = t[task_idx].high;
      if (i + n_valid < BISECT_K) {
        my_end = counts[i + n_valid];
        my_high = shifts[i + n_valid];
      }
      if (counts[i] <= my_end - 1) {
        if ((my_high - shifts[i]) / fabs(shifts[i]) > eps
            && my_high - shifts[i] > std::numeric_limits<double>::min()) {
          task_queue.push({counts[i], my_end, shifts[i], my_high});
        } else {
          int my_start = t[task_idx].start;
          if (i - n_valid >= 0) {
            my_start = counts[i - n_valid];
          }
          const int n_eq = my_end - counts[i];
          low.segment(counts[i], n_eq)
              = Eigen::VectorXd::Constant(n_eq, shifts[i]);
          high.segment(counts[i], n_eq)
              = Eigen::VectorXd::Constant(n_eq, my_high);
        }
      }
    }
  }
  low = low.reverse().eval();
  high = high.reverse().eval();
}

/**
 * Calculates an eigenvector from twisted factorization T - shift * I = L+ * D+
 * * L+^T = U- * D- * U-^T.
 * @param l_plus Subdiagonal of the L+.
 * @param u_minus Superdiagonal of the U-.
 * @param subdiagonal Subdiagonal of T
 * @param i At which column of `eigenvecs` to store resulting vector.
 * @param twist_idx Twist index.
 * @param[out] eigenvectors Matrix in which to store resulting vector.
 */
void calculate_eigenvector(const Eigen::VectorXd& l_plus,
                           const Eigen::VectorXd& u_minus,
                           const Eigen::Ref<const Eigen::VectorXd>& subdiagonal,
                           const int i, const int twist_idx,
                           Eigen::Ref<Eigen::MatrixXd>& eigenvectors) {
  using std::isinf;
  using std::isnan;
  auto vec = eigenvectors.col(i);
  const int n = vec.size();
  vec[twist_idx] = 1;
  for (int j = twist_idx + 1; j < n; j++) {
    if (vec[j - 1] != 0) {
      vec[j] = -u_minus[j - 1] * vec[j - 1];
    } else {
      vec[j] = -subdiagonal[j - 2] * vec[j - 2] / subdiagonal[j - 1];
      if (isnan(vec[j]) || isinf(vec[j])) {  // subdiagonal[j - 1]==0
        vec[j] = 0;
      }
    }
  }
  for (int j = twist_idx - 1; j >= 0; j--) {
    if (vec[j + 1] != 0) {
      vec[j] = -l_plus[j] * vec[j + 1];
    } else {
      vec[j] = -subdiagonal[j + 1] * vec[j + 2] / subdiagonal[j];
      if (isnan(vec[j]) || isinf(vec[j])) {  // subdiagonal[j]==0
        vec[j] = 0;
      }
    }
  }
  vec *= 1. / vec.norm();
}

/**
 * Finds good shift and shifts a LDL decomposition so as to keep element growth
 * low. L * D * L^T - shift * I = L2 * D2 * L2^T.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param low Low bound on wanted shift.
 * @param high High bound on wanted shift.
 * @param max_ele_growth Maximum desired element growth. If no better options
 * are found it might be exceeded.
 * @param max_shift Maximal difference of shhift from wanted bounds.
 * @param[out] l2 Subdiagonal of L2.
 * @param[out] d2 Diagonal of D2.
 * @param[out] shift Shift.
 * @param[out] min_element_growth Element growth achieved with resulting shift.
 */
void find_shift(const Eigen::VectorXd& l, const Eigen::VectorXd& d,
                const double low, const double high,
                const double max_ele_growth, const double max_shift,
                Eigen::VectorXd& l2, Eigen::VectorXd& d2, double& shift,
                double& min_element_growth) {
  Eigen::VectorXd l3(l2.size()), d3(d2.size());
  const std::vector<double> shifts = {
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
  for (double sh : shifts) {
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
double find_initial_shift(const Eigen::Ref<const Eigen::VectorXd> diagonal,
                          const Eigen::Ref<const Eigen::VectorXd> subdiagonal,
                          Eigen::VectorXd& l0, Eigen::VectorXd& d0,
                          const double min_eigval, const double max_eigval,
                          const double max_ele_growth) {
  double shift = 0;
  if (min_eigval > 0) {
    shift = 0.9 * min_eigval;
  } else if (max_eigval <= 0) {
    shift = 0.9 * max_eigval;
  }
  double element_growth = get_ldl(diagonal, subdiagonal, shift, l0, d0);
  if (element_growth < max_ele_growth) {
    return shift;
  }
  double plus = (max_eigval - min_eigval) * 1e-15;
  while (!(element_growth < max_ele_growth)) {  // if condition is flipped it
                                                // would be wrong for the case
                                                // where element_growth is nan
    plus *= -2;
    element_growth = get_ldl(diagonal, subdiagonal, shift + plus, l0, d0);
  }
  return shift + plus;
}

struct mrrr_task {
  int start, end;
  double shift;  // total shift, not just the last one
  Eigen::VectorXd l, d;
  int level;
};

/**
 * Calculates eigenvalues and eigenvectors of a irreducible tridiagonal matrix T
 * using MRRR algorithm. Use `tridiagonal_eigensolver` if any subdiagonal
 * element might be (very close to) zero.
 * @param diagonal Diagonal of of T.
 * @param subdiagonal Subdiagonal of T.
 * @param eigenvalues[out] Eigenvlues.
 * @param eigenvectors[out] Eigenvectors.
 * @param min_rel_sep Minimal relative separation of eigenvalues before
 * computing eigenvectors.
 * @param max_ele_growth Maximal desired element growth of LDL decompositions.
 */
void mrrr(const Eigen::Ref<const Eigen::VectorXd> diagonal,
          const Eigen::Ref<const Eigen::VectorXd> subdiagonal,
          Eigen::Ref<Eigen::VectorXd> eigenvalues,
          Eigen::Ref<Eigen::MatrixXd> eigenvectors,
          const double min_rel_sep = 1e-1, const double max_ele_growth = 2) {
  using std::copysign;
  using std::fabs;
  const double shift_error = 1e-14;
  const int n = diagonal.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  get_gresgorin(diagonal, subdiagonal, min_eigval, max_eigval);
  Eigen::VectorXd l(n - 1), d(n), l0(n - 1), d0(n);
  const double shift0 = find_initial_shift(
      diagonal, subdiagonal, l0, d0, min_eigval, max_eigval, max_ele_growth);
  for (int i = 0; i < n; i++) {
    if (i != n - 1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i] = d0[i] * get_random_perturbation_multiplier();
  }
  const Eigen::VectorXd subdiagonal_squared
      = subdiagonal.array() * subdiagonal.array();

  eigenvals_bisect(diagonal, subdiagonal_squared, min_eigval, max_eigval, low,
                   high);
  eigenvalues = (high + low) * 0.5;
  low.array() -= shift0;
  high.array() -= shift0;
  for (int i = 0; i < n; i++) {
    low[i] = low[i] * (1 - copysign(perturbation_range * n, low[i]));
    high[i] = high[i] * (1 + copysign(perturbation_range * n, high[i]));
    eigenval_bisect_refine(l, d, low[i], high[i], i);
  }
  std::queue<mrrr_task> block_queue;
  block_queue.push(mrrr_task{0, n, shift0, std::move(l), std::move(d), 0});
  l.resize(n - 1);  // after move out
  d.resize(n);
  while (!block_queue.empty()) {
    const mrrr_task block = block_queue.front();
    block_queue.pop();
    double shift = std::numeric_limits<double>::infinity();
    double min_element_growth = std::numeric_limits<double>::infinity();
    Eigen::VectorXd l2(n - 1), d2(n), l_plus(n - 1), u_minus(n - 1);
    for (int i = block.start; i < block.end; i++) {
      int cluster_end;
      for (cluster_end = i + 1; cluster_end < block.end; cluster_end++) {
        const int prev = cluster_end - 1;
        const double end_threshold
            = low[prev] * (1 - copysign(shift_error, low[prev]));
        if (high[cluster_end] < end_threshold) {
          break;
        }
      }
      cluster_end--;  // now this is the index of the last element of the
                      // cluster
      if (cluster_end - i > 0) {  // cluster
        const double max_shift = (high[i] - low[cluster_end]) * 10;
        double next_shift, min_ele_growth;
        find_shift(block.l, block.d, low[cluster_end], high[i], max_ele_growth,
                   max_shift, l, d, next_shift, min_ele_growth);
        for (int j = i; j <= cluster_end; j++) {
          low[j] = low[j] * (1 - copysign(shift_error, low[j])) - next_shift;
          high[j] = high[j] * (1 + copysign(shift_error, high[j])) - next_shift;
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
        const double low_gap = i == block.start
                                   ? std::numeric_limits<double>::infinity()
                                   : low[i - 1] - high[i];
        const double high_gap = i == block.end - 1
                                    ? std::numeric_limits<double>::infinity()
                                    : low[i] - high[i + 1];
        const double min_gap = std::min(low_gap, high_gap);
        const Eigen::VectorXd *l_ptr, *d_ptr;
        if (!(fabs(min_gap / ((high[i] + low[i]) * 0.5)) > min_rel_sep)) {
          if (!(fabs(min_gap / ((high[i] + low[i]) * 0.5 - shift)) > min_rel_sep
                && min_element_growth < max_ele_growth)) {
            const double max_shift = min_gap / min_rel_sep;
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
        twist_idx = get_twisted_factorization(
            *l_ptr, *d_ptr, (low[i] + high[i]) * 0.5, l_plus, u_minus);
        calculate_eigenvector(l_plus, u_minus, subdiagonal, i, twist_idx,
                              eigenvectors);
      }
    }
  }
}

/**
 * Calculates eigenvalues and eigenvectors of a tridiagonal matrix T using MRRR
 * algorithm. If a subdiagonal element is close to zero compared to neighbors on
 * diagonal the problem can be split into smaller ones.
 * @param diagonal Diagonal of of T.
 * @param subdiagonal Subdiagonal of T.
 * @param eigenvalues[out] Eigenvlues.
 * @param eigenvectors[out] Eigenvectors.
 * @param split_threshold Threshold for splitting the problem
 */
void tridiagonal_eigensolver(const Eigen::VectorXd& diagonal,
                             const Eigen::VectorXd& subdiagonal,
                             Eigen::VectorXd& eigenvalues,
                             Eigen::MatrixXd& eigenvectors,
                             const double split_threshold = 1e-12) {
  using std::fabs;
  const int n = diagonal.size();
  eigenvectors.resize(n, n);
  eigenvalues.resize(n);
  int last = 0;
  for (int i = 0; i < subdiagonal.size(); i++) {
    if (fabs(subdiagonal[i] / diagonal[i]) < split_threshold
        && fabs(subdiagonal[i] / diagonal[i + 1]) < split_threshold) {
      eigenvectors.block(last, i + 1, i + 1 - last, n - i - 1)
          = Eigen::MatrixXd::Constant(i + 1 - last, n - i - 1, 0);
      eigenvectors.block(i + 1, last, n - i - 1, i + 1 - last)
          = Eigen::MatrixXd::Constant(n - i - 1, i + 1 - last, 0);
      if (last == i) {
        eigenvectors(last, last) = 1;
        eigenvalues[last] = diagonal[last];
      } else {
        mrrr(diagonal.segment(last, i + 1 - last),
             subdiagonal.segment(last, i - last),
             eigenvalues.segment(last, i + 1 - last),
             eigenvectors.block(last, last, i + 1 - last, i + 1 - last));
      }

      last = i + 1;
    }
  }
  if (last == n - 1) {
    eigenvectors(last, last) = 1;
    eigenvalues[last] = diagonal[last];
  } else {
    mrrr(diagonal.segment(last, n - last),
         subdiagonal.segment(last, subdiagonal.size() - last),
         eigenvalues.segment(last, n - last),
         eigenvectors.block(last, last, n - last, n - last));
  }
}

}  // namespace internal
}  // namespace math
}  // namespace stan

#endif
