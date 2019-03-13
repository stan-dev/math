#ifndef STAN_MATH_PRIM_MAT_FUN_MRRR_HPP
#define STAN_MATH_PRIM_MAT_FUN_MRRR_HPP

#include <queue>
#include <cmath>

#include <Eigen/Dense>

namespace stan {
namespace math {

const double perturbation_range = 1e-15;

/**
 * Generates a random number for perturbing a relatively robust representation
 * @return A uniformly distributed random number between `1 - perturbation_range / 2` and `1 + perturbation_range / 2`.
 */
inline double get_random_perturbation_multiplier() {
  static const double rand_norm = perturbation_range / RAND_MAX;
  static const double almost_one = 1 - perturbation_range * 0.5;
  return almost_one + std::rand() * rand_norm;
}

/**
 * Calculates LDL decomposition of a shifted triagonal matrix T. D is diagonal, L is lower unit triangular (diagonal elements are 1,
 * all elements except diagonal and subdiagonal are 0),T - shift * I = L * D * L^T. Also calculates element growth of D: sum(std::fabs(D)) / std::fabs(sum(D)).
 * @param diag Diagonal of T
 * @param subdiag Subdiagonal of T.
 * @param shift Shift.
 * @param[out] l Subdiagonal of L.
 * @param[out] d_plus Diagonal of D.
 * @return Element growth.
 */
double get_ldl(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, double shift, Eigen::VectorXd& l, Eigen::VectorXd& d_plus) {
  d_plus[0] = diag[0] - shift;
  double element_growth = std::fabs(d_plus[0]);
  double element_growth_denominator = d_plus[0];
  for (int i = 0; i < subdiag.size(); i++) {
    l[i] = subdiag[i] / d_plus[i];
    d_plus[i] *= get_random_perturbation_multiplier();
    d_plus[i + 1] = diag[i + 1] - shift - l[i] * subdiag[i];
    l[i] *= get_random_perturbation_multiplier();
    element_growth += std::fabs(d_plus[i + 1]);
    element_growth_denominator += d_plus[i + 1];
  }
  d_plus[subdiag.size()] *= get_random_perturbation_multiplier();
  return element_growth / std::fabs(element_growth_denominator);
}

/**
 * Shifts a LDL decomposition. The algorithm is sometimes called stationary quotients-differences with shifts (stqds).
 * D and D+ are diagonal, L and L+ are lower unit triangular (diagonal elements are 1,
 * all elements except diagonal and subdiagonal are 0). L * D * L^T - shift * I = L+ * D * L+^T.
 * Also calculates element growth of D+: sum(std::fabs(D+)) / std::fabs(sum(D+)).
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param shift Shift.
 * @param[out] l_plus Subdiagonal of L+.
 * @param[out] d_plus Diagonal of D+.
 * @return Element growth.
 */
double get_shifted_ldl(const Eigen::VectorXd& l, const Eigen::VectorXd& d, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& d_plus) {
  int n = l.size();
  double s = -shift;
  double element_growth = 0;
  double element_growth_denominator = 0;
  for (int i = 0; i < n; i++) {
    d_plus[i] = s + d[i];
    element_growth += std::fabs(d_plus[i]);
    element_growth_denominator += d_plus[i];
    l_plus[i] = l[i] * (d[i] / d_plus[i]);
    if (std::isinf(d_plus[i]) && std::isinf(s)) { // this happens if d_plus[i]==0 -> in next iteration d_plus==inf and s==inf
      s = l[i] * l[i] * d[i] - shift;
    }
    else {
      s = l_plus[i] * l[i] * s - shift;
    }
  }
  d_plus[n] = s + d[n];
  element_growth += std::fabs(d_plus[n]);
  return element_growth / std::fabs(element_growth_denominator);
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
 * @return Twist index.
 */
int get_twisted_factorization(const Eigen::VectorXd& l, const Eigen::VectorXd& d, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& u_minus) {
  int n = l.size();
  //calculate shifted ldl
  Eigen::VectorXd s(n + 1);
  s[0] = -shift;
  for (int i = 0; i < n; i++) {
    double d_plus = s[i] + d[i];
    l_plus[i] = l[i] * (d[i] / d_plus);
    if (std::isnan(l_plus[i])) { //d_plus==0
      if (std::fabs(l[i]) < std::fabs(d[i])) { //one (or both) of d[i], l[i] is very close to 0
        l_plus[i] = d[i] * copysign(1., l[i]) * copysign(1., d_plus);
      }
      else {
        l_plus[i] = l[i] * copysign(1., d[i]) * copysign(1., d_plus);
      }
    }
    s[i + 1] = l_plus[i] * l[i] * s[i] - shift;
    if (std::isnan(s[i + 1])) {
      if (std::fabs(l_plus[i]) > std::fabs(s[i])) { //l_plus[i]==inf
        if (std::fabs(s[i]) > std::fabs(l[i])) { //l[i]==0
          s[i + 1] = s[i] * copysign(1., l[i]) * copysign(1., l_plus[i]) - shift;
        }
        else { //s[i]==0
          s[i + 1] = l[i] * copysign(1., s[i]) * copysign(1., l_plus[i]) - shift;
        }
      }
      else { //s[i]==inf
        if (std::fabs(l_plus[i]) > std::fabs(l[i])) { //l[i]==0
          s[i + 1] = l_plus[i] * copysign(1., l[i]) * copysign(1., s[i]) - shift;
        }
        else { //l_plus[i]==0
          s[i + 1] = l[i] * copysign(1., s[i]) * copysign(1., l_plus[i]) - shift;
        }
      }
    }
  }
  //calculate shifted udu and twist index
  double p = d[n] - shift;
  double min_gamma = std::fabs(s[n] + d[n]);
  int twist_index = n;

  for (int i = n - 1; i >= 0; i--) {
    double d_minus = d[i] * l[i] * l[i] + p;
    double t = d[i] / d_minus;
    u_minus[i] = l[i] * t;
    if (std::isnan(u_minus[i])) {
      if (std::isnan(t)) {
        t = copysign(1., d[i]) * copysign(1., d_minus);
        u_minus[i] = l[i] * t;
      }
      else { //t==inf, l[i]==0
        u_minus[i] = d[i] * copysign(1., l[i]) * copysign(1., t);
      }
    }
    double gamma = std::fabs(s[i] + t * p);
    if (std::isnan(gamma)) { //t==inf, p==0 OR t==0, p==inf
      double d_sign = d[i] * copysign(1., d_minus) * copysign(1., t);
      gamma = std::fabs(s[i] + d_sign);
      p = d_sign - shift;
    }
    else { //general case
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
 * Calculates Sturm count of a LDL decomposition of a tridiagonal matrix - number of eigenvalues larger or equal to shift.
 * Uses stqds - calculation of shifted LDL decomposition algorithm and counts number of positive elements in D.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param shift Shift.
 * @return Sturm count.
 */
int getSturmCountLdl(const Eigen::VectorXd& l, const Eigen::VectorXd& d, double shift) {
  int n = l.size();
  double s = -shift;
  double d_plus;
  int count = 0;
  for (int i = 0; i < n; i++) {
    d_plus = s + d[i];
    count += d_plus >= 0;
    if (std::isinf(d_plus) && std::isinf(s)) { // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
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
 * Refines bounds on the i-th largest eigenvalue of LDL decomposition of a matrix using bisection.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param low[in,out] Low bound on the eigenvalue.
 * @param high[in,out] High bound on the eigenvalue.
 * @param i i-th eigenvalue
 */
void eigenvalBisectRefine(const Eigen::VectorXd& l, const Eigen::VectorXd& d, double& low, double& high, int i) {
  double eps = 3e-16;
  while (std::fabs((high - low) / (high + low)) > eps && std::fabs(high - low) > std::numeric_limits<double>::min()) { // second term is for the case where the eigenvalue is 0 and division yields NaN
    double mid = (high + low) * 0.5;
    if (getSturmCountLdl(l, d, mid) > i) {
      low = mid;
    }
    else {
      high = mid;
    }
  }
}

/**
 * Calculates bounds on eigenvalues of a symmetric tridiagonal matrix T using Gresgorin discs.
 * @param diag Diagonal of T
 * @param subdiag Subdiagonal of T
 * @param min_eigval[out] Lower bound on eigenvalues.
 * @param max_eigval[out] Upper bound on eigenvalues.
 */
void getGresgorin(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, double& min_eigval, double& max_eigval) {
  int n = diag.size();
  min_eigval = diag[0] - std::fabs(subdiag[0]);
  max_eigval = diag[0] + std::fabs(subdiag[0]);
  for (int i = 1; i < n - 1; i++) {
    min_eigval = std::min(min_eigval, diag[i] - std::fabs(subdiag[i]) - std::fabs(subdiag[i - 1]));
    max_eigval = std::max(max_eigval, diag[i] + std::fabs(subdiag[i]) + std::fabs(subdiag[i - 1]));
  }
  min_eigval = std::min(min_eigval, diag[n - 1] - std::fabs(subdiag[n - 2]));
  max_eigval = std::max(max_eigval, diag[n - 1] + std::fabs(subdiag[n - 2]));
}

const int BISECT_K = 8;

/**
 * Calculates lower Sturm count of a tridiagonal matrix T - number of eigenvalues lower than shift for up to BISECT_K different shifts.
 * @param diag Diagonal of T.
 * @param subdiagSquared Squared elements of subdiagonal of T.
 * @param shifts Up to 8 different shifts. First `n_valid` are used.
 * @param n_valid How many Sturm counts to actually compute.
 * @return Array of Sturm counts of size BISECT_K. First `n_valid` are actual results.
 */
Eigen::Array<int, BISECT_K, 1> getSturmCountTVec2(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::VectorXd& subdiagSquared, Eigen::Array<double, BISECT_K, 1>& shifts, int n_valid) {
  Eigen::Array<double, BISECT_K, 1> d;
  d.head(n_valid) = diag[0] - shifts.head(n_valid);
  Eigen::Array<int, BISECT_K, 1> counts;
  counts.head(n_valid) = (d.head(n_valid) < 0).cast<int>();
  for (int j = 1; j < diag.size(); j++) {
    d.head(n_valid) = diag[j] - shifts.head(n_valid) - subdiagSquared[j - 1] / d.head(n_valid);
    counts.head(n_valid) += (d.head(n_valid) < 0).cast<int>();
  }
  return counts;
}

struct bisectionTask {
    int start, end;
    double low, high;
};

/**
 * Calculates eigenvalues of tridiagonal matrix T using bisection.
 * @param diag Diagonal of T.
 * @param subdiagSquared Squared elements of the subdiagonal.
 * @param min_eigval Lower bound on all eigenvalues.
 * @param max_eigval Upper bound on all eigenvalues.
 * @param low[out] Lower bounds on eigenvalues.
 * @param high[out] Upper bounds on eigenvalues.
 */
void eigenvalsBisect4(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::VectorXd& subdiagSquared, double min_eigval, double max_eigval, Eigen::VectorXd& low, Eigen::VectorXd& high) {
  int n = diag.size();
  double eps = 3e-16;

  std::queue<bisectionTask> tQueue;
  tQueue.push(bisectionTask{0, n, min_eigval, max_eigval});
  while (!tQueue.empty()) {
    int n_valid = std::min(BISECT_K, static_cast<int>(tQueue.size()));
    Eigen::Array<double, BISECT_K, 1> shifts;
    bisectionTask t[BISECT_K];
    for (int i = 0; i < n_valid; i++) {
      t[i] = tQueue.front();
      tQueue.pop();
    }
    for (int i = 0; i < BISECT_K; i++) {
      int task_idx = i % n_valid;
      int idx_in_task = i / n_valid;
      int task_total = BISECT_K / n_valid + (BISECT_K % n_valid > task_idx);
      shifts[i] = t[task_idx].low + (t[task_idx].high - t[task_idx].low) * (idx_in_task + 1.) / (task_total + 1);
    }
    Eigen::Array<int, BISECT_K, 1> counts = getSturmCountTVec2(diag, subdiagSquared, shifts, BISECT_K);
    for (int i = 0; i < n_valid; i++) {
      if (counts[i] >= t[i].start + 1) {
        if ((t[i].high - shifts[i]) / std::fabs(shifts[i]) > eps && shifts[i] - t[i].low > std::numeric_limits<double>::min()) {
          tQueue.push({t[i].start, counts[i], t[i].low, shifts[i]});
        }
        else {
          int n_eq = counts[i] - t[i].start;
          low.segment(t[i].start, n_eq) = Eigen::VectorXd::Constant(n_eq, t[i].low);
          high.segment(t[i].start, n_eq) = Eigen::VectorXd::Constant(n_eq, shifts[i]);
        }
      }
    }
    for (int i = 0; i < BISECT_K; i++) {
      int task_idx = i % n_valid;
      int idx_in_task = i / n_valid;
      int task_total = BISECT_K / n_valid + (BISECT_K % n_valid > task_idx);
      int my_end = t[task_idx].end;
      double my_high = t[task_idx].high;
      if (i + n_valid < BISECT_K) {
        my_end = counts[i + n_valid];
        my_high = shifts[i + n_valid];
      }
      if (counts[i] <= my_end - 1) {
        if ((my_high - shifts[i]) / std::fabs(shifts[i]) > eps && my_high - shifts[i] > std::numeric_limits<double>::min()) {
          tQueue.push({counts[i], my_end, shifts[i], my_high});
        }
        else {
          int my_start = t[task_idx].start;
          if (i - n_valid >= 0) {
            my_start = counts[i - n_valid];
          }
          int n_eq = my_end - counts[i];
          low.segment(counts[i], n_eq) = Eigen::VectorXd::Constant(n_eq, shifts[i]);
          high.segment(counts[i], n_eq) = Eigen::VectorXd::Constant(n_eq, my_high);
        }
      }
    }
  }
  low = low.reverse().eval();
  high = high.reverse().eval();
}

/**
 * Calculates an eigenvector from twisted factorization T - shift * I = L+ * D+ * L+^T = U- * D- * U-^T.
 * @param l_plus Subdiagonal of the L+.
 * @param u_minus Superdiagonal of the U-.
 * @param subdiag Subdiagonal of T
 * @param i At which column of `eigenvecs` to store resulting vector.
 * @param twist_idx Twist index.
 * @param[out] eigenvecs Matrix in which to store resulting vector.
 */
void calculateEigenvector(const Eigen::VectorXd& l_plus, const Eigen::VectorXd& u_minus, const Eigen::Ref<const Eigen::VectorXd>& subdiag, int i, int twist_idx, Eigen::Ref<Eigen::MatrixXd>& eigenvecs) {
  auto vec = eigenvecs.col(i);
  int n = vec.size();
  vec[twist_idx] = 1;
  for (int j = twist_idx + 1; j < n; j++) {
    if (vec[j - 1] != 0) {
      vec[j] = -u_minus[j - 1] * vec[j - 1];
    }
    else {
      vec[j] = -subdiag[j - 2] * vec[j - 2] / subdiag[j - 1];
      if (std::isnan(vec[j]) || std::isinf(vec[j])) { //subdiag[j - 1]==0
        vec[j] = 0;
      }
    }
  }
  for (int j = twist_idx - 1; j >= 0; j--) {
    if (vec[j + 1] != 0) {
      vec[j] = -l_plus[j] * vec[j + 1];
    }
    else {
      vec[j] = -subdiag[j + 1] * vec[j + 2] / subdiag[j];
      if (std::isnan(vec[j]) || std::isinf(vec[j])) { //subdiag[j]==0
        vec[j] = 0;
      }
    }
  }
  vec *= 1./vec.norm();
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
 * @param[out] shift Shift.
 * @param[out] min_element_growth Element growth achieved with resulting shift.
 */
void findShift(const Eigen::VectorXd& l, const Eigen::VectorXd& d, double low, double high, double max_ele_growth, double max_shift,
               Eigen::VectorXd& l2, Eigen::VectorXd& d2, double& shift, double& min_element_growth) {
  Eigen::VectorXd l3(l2.size()), d3(d2.size());
  std::vector<double> shifts = {
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
    double element_growth = get_shifted_ldl(l, d, sh, l3, d3);
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
 * Finds a good value for shift of the initial LDL factorization T - shift * I = L * D * L^T.
 * @param diag Diagonal of T.
 * @param subdiag Subdiagonal of T.
 * @param l0 Subdiagonal of L.
 * @param d0 Diagonal of D.
 * @param min_eigval Lower bound on eigenvalues of T.
 * @param max_eigval High bound on eigenvalues of T
 * @param max_ele_growth Maximum desired element growth.
 * @return
 */
double findInitialShift(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, Eigen::VectorXd& l0, Eigen::VectorXd& d0, double min_eigval, double max_eigval,
                        double max_ele_growth) {
  double shift = 0;
  if (min_eigval > 0) {
    shift = 0.9 * min_eigval;
  }
  else if (max_eigval <= 0) {
    shift = 0.9 * max_eigval;
  }
  double element_growth = get_ldl(diag, subdiag, shift, l0, d0);
  if (element_growth < max_ele_growth) {
    return shift;
  }
  double plus = (max_eigval - min_eigval) * 1e-15;
  while (!(element_growth < max_ele_growth)) { //if condition is flipped it would be wrong for the case where element_growth is nan
    plus *= -2;
    element_growth = get_ldl(diag, subdiag, shift + plus, l0, d0);
  }
  return shift + plus;
}

struct mrrrTask {
    int start, end;
    double shift; //total shift, not just last one
    Eigen::VectorXd l, d;
    int level;
};

/**
 * Calculates eigenvalues and eigenvectors of a irreducible tridiagonal matrix T using MRRR algorithm. Use `tridiagonal_eigensolver` if any subdiagonal element bight be (very close to) zero.
 * @param diag Diagonal of of T.
 * @param subdiag Subdiagonal of T.
 * @param eigenvals[out] Eigenvlues.
 * @param eigenvecs[out] Eigenvectors.
 * @param min_rel_sep Minimal relative separation of eigenvalues before computing eigenvectors.
 * @param max_ele_growth Maximal desired element growth of LDL decompositions.
 */
void mrrr(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, Eigen::Ref<Eigen::VectorXd> eigenvals,
           Eigen::Ref<Eigen::MatrixXd> eigenvecs, double min_rel_sep = 1e-1, double max_ele_growth = 2) {
  double shift_error = 1e-14;
  int n = diag.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  getGresgorin(diag, subdiag, min_eigval, max_eigval);
  Eigen::VectorXd l(n - 1), d(n), l0(n - 1), d0(n);
  double shift0 = findInitialShift(diag, subdiag, l0, d0, min_eigval, max_eigval, max_ele_growth);
  for (int i = 0; i < n; i++) {
    if (i != n - 1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i] = d0[i] * get_random_perturbation_multiplier();
  }
  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();

  eigenvalsBisect4(diag, subdiagSquared, min_eigval, max_eigval, low, high);
  eigenvals = (high + low) * 0.5;
  low.array() -= shift0;
  high.array() -= shift0;
  for (int i = 0; i < n; i++) {
    low[i] = low[i] * (1 - copysign(perturbation_range * n, low[i]));
    high[i] = high[i] * (1 + copysign(perturbation_range * n, high[i]));
    eigenvalBisectRefine(l, d, low[i], high[i], i);
  }
  std::queue<mrrrTask> blockQueue;
  blockQueue.push(mrrrTask{0, n, shift0, std::move(l), std::move(d), 0});
  l.resize(n - 1); //after move out
  d.resize(n);
  while (!blockQueue.empty()) {
    mrrrTask block = blockQueue.front();
    blockQueue.pop();
    double shift = std::numeric_limits<double>::infinity();
    double min_element_growth = std::numeric_limits<double>::infinity();
    Eigen::VectorXd l2(n - 1), d2(n), l_plus(n - 1), u_minus(n - 1);
    for (int i = block.start; i < block.end; i++) {
      int cluster_end;
      for (cluster_end = i + 1; cluster_end < block.end; cluster_end++) {
        int prev = cluster_end - 1;
        double end_threshold = low[prev] * (1 - copysign(shift_error, low[prev]));
        if (high[cluster_end] < end_threshold) {
          break;
        }
      }
      cluster_end--; //now this is the index of the last element of the cluster
      if (cluster_end - i > 0) {//cluster
        double max_shift = (high[i] - low[cluster_end]) * 10;
        double currentShift, min_ele_growth;
        findShift(block.l, block.d, low[cluster_end], high[i], max_ele_growth, max_shift, l, d, currentShift, min_ele_growth);
        for (int j = i; j <= cluster_end; j++) {
          low[j] = low[j] * (1 - copysign(shift_error, low[j])) - currentShift;
          high[j] = high[j] * (1 + copysign(shift_error, high[j])) - currentShift;
          eigenvalBisectRefine(l, d, low[j], high[j], j);
        }
        blockQueue.push(mrrrTask{i, cluster_end + 1, block.shift + currentShift, std::move(l), std::move(d), block.level + 1});
        l.resize(n - 1); //after move out
        d.resize(n);

        i = cluster_end;
      }
      else { //isolated eigenvalue
        double low_t = low[i], high_t = high[i];
        int twist_idx;
        double low_gap = i == block.start ? std::numeric_limits<double>::infinity() : low[i - 1] - high[i];
        double high_gap = i == block.end - 1 ? std::numeric_limits<double>::infinity() : low[i] - high[i + 1];
        double min_gap = std::min(low_gap, high_gap);
        Eigen::VectorXd *l_ptr, *d_ptr;
        if(!(std::fabs(min_gap / ((high[i] + low[i]) * 0.5)) > min_rel_sep)){
          if (!(std::fabs(min_gap / ((high[i] + low[i]) * 0.5 - shift)) > min_rel_sep && min_element_growth < max_ele_growth)){
            double max_shift = min_gap / min_rel_sep;
            findShift(block.l, block.d, low[i], high[i], max_ele_growth, max_shift, l2, d2, shift, min_element_growth);
          }
          low[i] = low[i] * (1 - copysign(shift_error, low[i])) - shift;
          high[i] = high[i] * (1 + copysign(shift_error, high[i])) - shift;
          eigenvalBisectRefine(l2, d2, low[i], high[i], i);
          l_ptr=&l2;
          d_ptr=&d2;
        }
        else{
          l_ptr=&block.l;
          d_ptr=&block.d;
        }
        twist_idx = get_twisted_factorization(*l_ptr, *d_ptr, (low[i] + high[i]) * 0.5, l_plus, u_minus);
        calculateEigenvector(l_plus, u_minus, subdiag, i, twist_idx, eigenvecs);
      }
    }
  }
}

/**
* Calculates eigenvalues and eigenvectors of a tridiagonal matrix T using MRRR algorithm.
* If a subdiagonal element is close to zero compared to neighbors on diagonal the problem can be split into smaller ones.
* @param diag Diagonal of of T.
* @param subdiag Subdiagonal of T.
* @param eigenvals[out] Eigenvlues.
* @param eigenvecs[out] Eigenvectors.
* @param splitThreshold Threshold for splitting the problem
*/
void tridiagonal_eigensolver(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs, double splitThreshold = 1e-12) {
  int n = diag.size();
  eigenvecs.resize(n, n);
  eigenvals.resize(n);
  int last = 0;
  for (int i = 0; i < subdiag.size(); i++) {
    if (std::fabs(subdiag[i] / diag[i]) < splitThreshold && std::fabs(subdiag[i] / diag[i + 1]) < splitThreshold) {
      eigenvecs.block(last, i + 1, i + 1 - last, n - i - 1) = Eigen::MatrixXd::Constant(i + 1 - last, n - i - 1, 0);
      eigenvecs.block(i + 1, last, n - i - 1, i + 1 - last) = Eigen::MatrixXd::Constant(n - i - 1, i + 1 - last, 0);
      if (last == i) {
        eigenvecs(last, last) = 1;
        eigenvals[last] = diag[last];
      }
      else {
        mrrr(diag.segment(last, i + 1 - last),
                 subdiag.segment(last, i - last),
                 eigenvals.segment(last, i + 1 - last),
                 eigenvecs.block(last, last, i + 1 - last, i + 1 - last));
      }

      last = i + 1;
    }
  }
  if (last == n - 1) {
    eigenvecs(last, last) = 1;
    eigenvals[last] = diag[last];
  }
  else {
    mrrr(diag.segment(last, n - last),
             subdiag.segment(last, subdiag.size() - last),
             eigenvals.segment(last, n - last),
             eigenvecs.block(last, last, n - last, n - last));
  }
}

}
}

#endif
