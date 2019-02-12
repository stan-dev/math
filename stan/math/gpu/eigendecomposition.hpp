//
// Created by tadej on 7. 12. 2018.
//

#include <Eigen/QR>
#include <iostream>
#include <queue>

#include <stan/math/gpu/matrix_gpu.hpp>

using namespace std;

#ifndef STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP
#define STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP

#ifdef STAN_OPENCL

//#define TIME_IT
//#define SKIP_Q

namespace stan {
namespace math {

void s(const Eigen::MatrixXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << ")" << std::endl;
}

void s(const Eigen::VectorXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << ")" << std::endl;
}

void s(const Eigen::RowVectorXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << ")" << std::endl;
}

void p(const Eigen::MatrixXd& a) {
  s(a);
  std::cout << a << std::endl;
}


void p(const Eigen::VectorXd& a) {
  s(a);
  std::cout << a << std::endl;
}

void p(const Eigen::RowVectorXd& a) {
  s(a);
  std::cout << a << std::endl;
}

void p(const matrix_gpu& a) {
  Eigen::MatrixXd b(a.rows(), a.cols());
  copy(b, a);
  s(b);
  std::cout << b << std::endl;
}

void block_householder_tridiag3(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
  packed = A;
  for (size_t k = 0; k < packed.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd V(packed.rows() - k - 1, actual_r);
    Eigen::MatrixXd U(packed.rows() - k - 1, actual_r);
    V.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);
    U.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(U.rows(), U.cols(), 0);

    for (size_t j = 0; j < actual_r; j++) {
      auto householder = packed.col(k + j).tail(packed.rows() - k - j - 1);
      if (j != 0) {
        auto householder_whole = packed.col(k + j).tail(packed.rows() - k - j);
        householder_whole -= U.block(j - 1, 0, householder_whole.size(), j) * V.block(j - 1, 0, 1, j).transpose() +
                             V.block(j - 1, 0, householder_whole.size(), j) * U.block(j - 1, 0, 1, j).transpose();
      }
      double q = householder.squaredNorm();
      double alpha = -copysign(sqrt(q), packed(k + j, k + j));
      q -= householder[0] * householder[0];
      householder[0] -= alpha;
      q += householder[0] * householder[0];
      q = sqrt(q);
      householder *= SQRT_2 / q;

      auto& u = householder;
      Eigen::VectorXd v(householder.size() + 1);
      v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - j - 1, packed.cols() - k - j - 1).selfadjointView<Eigen::Lower>() * u
                                   - U.bottomLeftCorner(u.size(), j) * (V.bottomLeftCorner(u.size(), j).transpose() * u)
                                   - V.bottomLeftCorner(u.size(), j) * (U.bottomLeftCorner(u.size(), j).transpose() * u);
      v[0] = q / SQRT_2;// - alpha * householder[1];
      double cnst = v.tail(householder.size()).transpose() * u;
      v.tail(householder.size()) -= 0.5 * cnst * u;

      packed(k + j, k + j + 1) = packed(k + j + 1, k + j) * q / SQRT_2 + alpha - v[0] * u[0];
      U.col(j).tail(U.rows() - j) = u;
      V.col(j).tail(V.rows() - j) = v.tail(V.rows() - j);
    }

    Eigen::MatrixXd& Y = U;
    Eigen::MatrixXd Y_partial_update = U * V.transpose();
    packed.block(k + actual_r, k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -=
            (Y_partial_update + Y_partial_update.transpose()).bottomRightCorner(Y_partial_update.rows() - actual_r + 1, Y_partial_update.cols() - actual_r + 1);
  }
  packed(packed.rows() - 2, packed.cols() - 1) = packed(packed.rows() - 1, packed.cols() - 2);
}

void block_apply_packed_Q3(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100) {
  //if input A==Identity, constructs Q
  Eigen::MatrixXd scratchSpace(A.rows(), r);
  for (int k = (packed.rows() - 3) / r * r; k >= 0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k - 1, actual_r);
    W.col(0) = packed.col(k).tail(W.rows());
    for (size_t j = 1; j < actual_r; j++) {
      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
      W.col(j).noalias() = -W.leftCols(j) * scratchSpace.col(0).head(j);
      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
    }
    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
  }
}

const double perturbation_range = 1e-14;

inline double get_random_perturbation_multiplier() {
  static const double rand_norm = perturbation_range / RAND_MAX;
  static const double almost_one = 1 - perturbation_range * 0.5;
  return almost_one + std::rand() * rand_norm;
}

void get_ldl(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, double shift, Eigen::VectorXd& l, Eigen::VectorXd& d_plus) {
  d_plus[0] = diag[0] - shift;
  for (int i = 0; i < subdiag.size(); i++) {
    l[i] = subdiag[i] / d_plus[i];
    //d_plus[i] *= get_random_perturbation_multiplier();
    d_plus[i + 1] = diag[i + 1] - shift - l[i] * subdiag[i];
    //l[i] *= get_random_perturbation_multiplier();
  }
  //d_plus[subdiag.size()] *= get_random_perturbation_multiplier();
}

//stationary qds
double get_perturbed_shifted_ldl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& d_plus) {
  int n = l.size();
  double s = -shift;
  double element_growth = 0;
  double element_growth_denominator = 0;
  for (int i = 0; i < n; i++) {
    double di_perturbed = d[i];// * get_random_perturbation_multiplier();
    double li_perturbed = l[i];// * get_random_perturbation_multiplier();
    d_plus[i] = s + di_perturbed;
    element_growth += abs(d_plus[i]);
    element_growth_denominator += d_plus[i];
    l_plus[i] = li_perturbed * (di_perturbed / d_plus[i]);
    if (is_inf(d_plus[i]) && is_inf(s)) { // this happens if d_plus[i]==0 -> in next iteration d_plus==inf and s==inf
      s = li_perturbed * li_perturbed * di_perturbed - shift;
    }
    else {
      s = l_plus[i] * li_perturbed * s - shift;
    }
  }
  d_plus[n] = s + d[n];// * get_random_perturbation_multiplier();
  element_growth += abs(d_plus[n]);
  return element_growth / abs(element_growth_denominator);
}

//dtwqds
int get_twisted_factorization(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& u_minus) {
  int n = l.size();
  //calculate shifted ldl
  Eigen::VectorXd s(n + 1);
  s[0] = -shift;
  for (int i = 0; i < n; i++) {
    double d_plus = s[i] + d[i];
    l_plus[i] = l[i] * (d[i] / d_plus);
    if (is_nan(l_plus[i])) { //d_plus==0
      //one (or both) of d[i], l[i] is very close to 0
      if (abs(l[i]) < abs(d[i])) {
        l_plus[i] = d[i] * copysign(1., l[i]) * copysign(1., d_plus);
      }
      else {
        l_plus[i] = l[i] * copysign(1., d[i]) * copysign(1., d_plus);
      }
    }
    s[i + 1] = l_plus[i] * l[i] * s[i] - shift;
    if (is_nan(s[i + 1])) {
      if (abs(l_plus[i]) > abs(s[i])) { //l_plus[i]==inf
        if (abs(s[i]) > abs(l[i])) { //l[i]==0
          s[i + 1] = s[i] * copysign(1., l[i]) * copysign(1., l_plus[i]) - shift;
        }
        else { //s[i]==0
          s[i + 1] = l[i] * copysign(1., s[i]) * copysign(1., l_plus[i]) - shift;
        }
      }
      else { //s[i]==inf
        if (abs(l_plus[i]) > abs(l[i])) { //l[i]==0
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
  double min_gamma = abs(s[n] + d[n]);
  int twist_index = n;

  for (int i = n - 1; i >= 0; i--) {
    double d_minus = d[i] * l[i] * l[i] + p;
    double t = d[i] / d_minus;
    u_minus[i] = l[i] * t;
    if (is_nan(u_minus[i])) {
      if (is_nan(t)) {
        t = copysign(1., d[i]) * copysign(1., d_minus);
        u_minus[i] = l[i] * t;
      }
      else { //t==inf, l[i]==0
        u_minus[i] = d[i] * copysign(1., l[i]) * copysign(1., t);
      }
    }
    double gamma = abs(s[i] + t * p);//TODO: all inf/nan -> need other shift!
    if (is_nan(gamma)) { //t==inf, p==0 OR t==0, p==inf
      double d_sign = d[i] * copysign(1., d_minus) * copysign(1., t);
      gamma = abs(s[i] + d_sign);
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

int getSturmCountLdlLdl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift) {
  int n = l.size();
  double s = -shift;
  double l_plus, d_plus;
  int count = 0;
  for (int i = 0; i < n; i++) {
    d_plus = s + d[i];
    count += d_plus >= 0;
    if (is_inf(d_plus) && is_inf(s)) { // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
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

double eigenvalBisectRefine(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double& low, double& high, int i) {
  int n = d.size();
  double eps = 3e-16;
  while (abs((high - low) / (high + low)) > eps && abs(high - low) > std::numeric_limits<double>::min()) { // second term is for the case where the eigenvalue is 0 and division yields NaN
    double mid = (high + low) * 0.5;
    if (getSturmCountLdlLdl(d, l, mid) > i) {
      low = mid;
    }
    else {
      high = mid;
    }
  }
}

void getGresgorin(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, double& min_eigval, double& max_eigval) {
  int n = diag.size();
  min_eigval = diag[0] - abs(subdiag[0]);
  max_eigval = diag[0] + abs(subdiag[0]);
  for (int i = 1; i < n - 1; i++) {
    min_eigval = std::min(min_eigval, diag[i] - abs(subdiag[i]) - abs(subdiag[i - 1]));
    max_eigval = std::max(max_eigval, diag[i] + abs(subdiag[i]) + abs(subdiag[i - 1]));
  }
  min_eigval = std::min(min_eigval, diag[n - 1] - abs(subdiag[n - 2]));
  max_eigval = std::max(max_eigval, diag[n - 1] + abs(subdiag[n - 2]));
}

const int BISECT_K = 4;

Eigen::Array<int, BISECT_K, 1> getSturmCountTVec(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiagSquared, Eigen::Array<double, BISECT_K, 1>& shifts, int n_valid) {
  Eigen::Array<double, BISECT_K, 1> d;
  d.head(n_valid) = diag[0] - shifts.head(n_valid);
  Eigen::Array<int, BISECT_K, 1> counts;
  counts.head(n_valid) = (d.head(n_valid) >= 0).cast<int>();
  for (int j = 1; j < diag.size(); j++) { //Sturm count via LDL factorization
    d.head(n_valid) = diag[j] - shifts.head(n_valid) - subdiagSquared[j - 1] / d.head(n_valid);
    counts.head(n_valid) += (d.head(n_valid) >= 0).cast<int>();
  }
  return counts;
}

Eigen::Array<int, BISECT_K, 1> getSturmCountTVec2(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiagSquared, Eigen::Array<double, BISECT_K, 1>& shifts, int n_valid) {
  Eigen::Array<double, BISECT_K, 1> d;
  d.head(n_valid) = diag[0] - shifts.head(n_valid);
  Eigen::Array<int, BISECT_K, 1> counts;
  counts.head(n_valid) = (d.head(n_valid) < 0).cast<int>();
  for (int j = 1; j < diag.size(); j++) { //Sturm count via LDL factorization
    d.head(n_valid) = diag[j] - shifts.head(n_valid) - subdiagSquared[j - 1] / d.head(n_valid);
    counts.head(n_valid) += (d.head(n_valid) < 0).cast<int>();
  }
  return counts;
}

void eigenvalsBisect3(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiagSquared, double min_eigval, double max_eigval, Eigen::VectorXd& low, Eigen::VectorXd& high) {
  int n = diag.size();
  double eps = 3e-16;

  for (int i = 0; i < n; i += BISECT_K) {
    int n_valid = std::min(BISECT_K, n - i);
    auto low_work = low.segment(i, n_valid).array();
    auto high_work = high.segment(i, n_valid).array();
    low_work = Eigen::Array<double, Eigen::Dynamic, 1>::Constant(n_valid, min_eigval);
    high_work = Eigen::Array<double, Eigen::Dynamic, 1>::Constant(n_valid, max_eigval);
    while ((abs((high_work - low_work) / low_work) > eps && abs(high_work - low_work) > std::numeric_limits<double>::min()).any()) {
      Eigen::Array<double, BISECT_K, 1> shifts;
      shifts.head(n_valid) = (high_work + low_work) * 0.5;
      Eigen::Array<int, BISECT_K, 1> counts = getSturmCountTVec(diag, subdiagSquared, shifts, n_valid);
      for (int k = 0; k < n_valid; k++) {
        if (counts[k] > i + k) {
          low_work[k] = shifts[k];
        }
        else {
          high_work[k] = shifts[k];
        }
      }
    }
  }
}

struct task{
    int start, end;
    double low, high;
};

void eigenvalsBisect4(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiagSquared, double min_eigval, double max_eigval, Eigen::VectorXd& low, Eigen::VectorXd& high) {
  int n = diag.size();
  double eps = 3e-16;

  std::queue<task> tQueue;
  tQueue.push(task{0, n, min_eigval, max_eigval});
  while(!tQueue.empty()){
    int n_valid = std::min(BISECT_K, static_cast<int>(tQueue.size()));
    Eigen::Array<double, BISECT_K, 1> shifts;
    task t[BISECT_K];
    for(int i=0;i<n_valid;i++){
      t[i] = tQueue.front();
      tQueue.pop();
    }
    for(int i=0;i<BISECT_K;i++) {
      int task_idx = i%n_valid;
      int idx_in_task = i/n_valid;
      int task_total = BISECT_K/n_valid + (BISECT_K%n_valid > task_idx);
      shifts[i] = t[task_idx].low + (t[task_idx].high - t[task_idx].low) * (idx_in_task+1.) / (task_total+1);
    }
    Eigen::Array<int, BISECT_K, 1> counts = getSturmCountTVec2(diag, subdiagSquared, shifts, BISECT_K);
    for(int i=0;i<n_valid;i++){
      if(counts[i]>=t[i].start+1){ //TODO +/- ?
        if((t[i].high - shifts[i])/abs(shifts[i]) > eps && abs(t[i].high - shifts[i]) > std::numeric_limits<double>::min()){
          tQueue.push({t[i].start, counts[i], t[i].low, shifts[i]});
        }
        else{
          int n_eq = counts[i] - t[i].start;
          low.segment(t[i].start, n_eq) = Eigen::VectorXd::Constant(n_eq, t[i].low);
          high.segment(t[i].start, n_eq) = Eigen::VectorXd::Constant(n_eq, shifts[i]);
        }
      }
    }
    for(int i=0;i<BISECT_K;i++) {
      int task_idx = i%n_valid;
      int idx_in_task = i/n_valid;
      int task_total = BISECT_K/n_valid + (BISECT_K%n_valid > task_idx);
      int my_end = t[task_idx].end;
      double my_high = t[task_idx].high;
      if(i+n_valid<BISECT_K){
        my_end = counts[i+n_valid];
        my_high = shifts[i+n_valid];
      }
      if(counts[i] <= my_end - 1){ //TODO +/- ?
        if((my_high - shifts[i])/abs(shifts[i]) > eps && abs(my_high - shifts[i]) > std::numeric_limits<double>::min()) {
          tQueue.push({counts[i], my_end, shifts[i], my_high});
        }
        else{
          int my_start = t[task_idx].start;
          if(i-n_valid>=0){
            my_start = counts[i-n_valid];
          }
          int n_eq = my_end - counts[i];
          low.segment(counts[i], n_eq) = Eigen::VectorXd::Constant(n_eq, shifts[i]);
          high.segment(counts[i], n_eq) = Eigen::VectorXd::Constant(n_eq, my_high);
        }
      }
    }
  }
  low=low.reverse().eval();
  high=high.reverse().eval();
}


void constructAndRefineCluster(const Eigen::VectorXd& l, const Eigen::VectorXd& d, int i, Eigen::VectorXd& high, Eigen::VectorXd& low, int& cluster_start, int& cluster_end) {
  int n = d.size();
  cluster_start = i;
  for (cluster_end = i; cluster_end < n - 1; cluster_end++) {
    int next = cluster_end + 1;
    low[next] = low[next] * (1 - copysign(perturbation_range * n, low[next]));
    high[next] = high[next] * (1 + copysign(perturbation_range * n, high[next]));
    eigenvalBisectRefine(d, l, low[next], high[next], next);
    double end_threshold = low[cluster_end] * (1 - copysign(sqrt(perturbation_range) * n, low[cluster_end]));
    if (high[next] < end_threshold) {
      break;
    }
  }
}


bool isClusterSeparable(const Eigen::VectorXd& high, const Eigen::VectorXd& low, int cluster_start, int cluster_end) {
  for (int j = cluster_start; j < cluster_end; j++) {
    if (high[j + 1] >= low[j]) {//problem - need new perturb
      return false;
    }
  }
  return true;
}

void perturbRepresentation(const Eigen::VectorXd& l0, const Eigen::VectorXd& d0, Eigen::VectorXd& l, Eigen::VectorXd& d) {
  int n = l0.size();
  for (int j = 0; j < n; j++) {
    l[j] = l0[j] * get_random_perturbation_multiplier();
    d[j] = d0[j] * get_random_perturbation_multiplier();
  }
  d[n] = d0[n] * get_random_perturbation_multiplier();
}

void calculateEigenvector(const Eigen::VectorXd& l_plus, const Eigen::VectorXd& u_minus, const Eigen::VectorXd& subdiag, int i, int twist_idx, Eigen::MatrixXd& eigenvecs) {
  auto vec = eigenvecs.col(i);
  int n = vec.size();
  vec[twist_idx] = 1;
  for (int j = twist_idx + 1; j < n; j++) {
    if (vec[j - 1] != 0) {
      vec[j] = -u_minus[j - 1] * vec[j - 1];
    }
    else {
      vec[j] = -subdiag[j - 2] * vec[j - 2] / subdiag[j - 1];
      if (is_nan(vec[j]) || is_inf(vec[j])) { //subdiag[j - 1]==0
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
      if (is_nan(vec[j]) || is_inf(vec[j])) { //subdiag[j]==0
        vec[j] = 0;
      }
    }
  }
  //vec /= vec.norm();
}

void findShift(const Eigen::VectorXd& l, const Eigen::VectorXd& d, double low, double high, double max_ele_growth, double max_shift,
        Eigen::VectorXd& l2, Eigen::VectorXd& d2, double& shift, double& min_element_growth) {
  Eigen::VectorXd l3(l2.size()), d3(d2.size());
  vector<double> shifts = {
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
  min_element_growth = numeric_limits<double>::infinity();
  for (double sh : shifts) {
    //sh -= shift0;
    double element_growth = get_perturbed_shifted_ldl(d, l, sh, l3, d3);
    cout << " element growth: " << element_growth << " at " << sh << endl;
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
  cout << "\t\t" << " element growth: " << min_element_growth << " at " << shift << endl;
}

void mrrr(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs, double min_rel_sep = 1e-2, double max_ele_growth = 2) {
  int n = diag.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  getGresgorin(diag, subdiag, min_eigval, max_eigval);
  //double span = max_eigval - min_eigval;
  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();
  auto start = std::chrono::steady_clock::now();
  eigenvalsBisect4(diag, subdiagSquared, min_eigval, max_eigval, low, high);
  cout << "bisect: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  Eigen::VectorXd l(n - 1), d(n), l0(n - 1), d0(n);
  double shift0 = min_eigval - (max_eigval - min_eigval) * 0.1;
  get_ldl(diag, subdiag, shift0, l0, d0);
  Eigen::VectorXd l2(n - 1), d2(n), l3(n - 1), d3(n), l_plus(n - 1), d_plus(n), u_minus(n - 1), d_minus(n), s(n), p(n);

  for (int i = 0; i < n; i++) {
    eigenvals[i] = (high[i] + low[i]) * 0.5;
    //from here on eigenvals[i] constins original eigenval
    if (i != n - 1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i] = d0[i] * get_random_perturbation_multiplier();
    low[i] -= shift0;
    high[i] -= shift0;
    //from here on high[i] and low[i] contain bounds on shifted representation
  }
  low[0] = low[0] * (1 - copysign(perturbation_range * n, low[0]));
  high[0] = high[0] * (1 + copysign(perturbation_range * n, high[0]));
  eigenvalBisectRefine(d, l, low[0], high[0], 0);

  double shift = std::numeric_limits<double>::infinity();
  double min_element_growth = std::numeric_limits<double>::infinity();
  int cluster_start = -1, cluster_end = -1;
  for (int i = 0; i < n; i++) {
    /*if (i > cluster_end) {
      constructAndRefineCluster(l, d, i, high, low, cluster_start, cluster_end);
      while (true) {
        // perturb the representation if any eigenvalues in the cluster are too close
        if (isClusterSeparable(high, low, cluster_start, cluster_end)) {
          break;
        }
        cout << "perturbing at i = " << i << " cluster: " << cluster_start << ", " << cluster_end << endl;
        perturbRepresentation(l0, d0, l, d);
        for (int j = cluster_start; j <= cluster_end; j++) {
          low[j] = low[j] * (1 - copysign(perturbation_range * n, low[j]));
          high[j] = high[j] * (1 + copysign(perturbation_range * n, high[j]));
          eigenvalBisectRefine(d, l, low[j], high[j], j);
        }
      }
      if (cluster_end != n - 1) {
        int next = cluster_end + 1;
        low[next] = low[next] * (1 - copysign(perturbation_range * n, low[next]));
        high[next] = high[next] * (1 + copysign(perturbation_range * n, high[next]));
        eigenvalBisectRefine(d, l, low[next], high[next], next);
      }
      min_element_growth = std::numeric_limits<double>::infinity();
    }*/
    if(i!=n-1){
      low[i] = low[i] * (1 - copysign(perturbation_range * n, low[i]));
      high[i] = high[i] * (1 + copysign(perturbation_range * n, high[i]));
      eigenvalBisectRefine(d, l, low[i+1], high[i+1], i+1);
    }
    int twist_idx;
    double low_gap = i == 0 ? std::numeric_limits<double>::infinity() : low[i - 1] - high[i];
    double high_gap = i == n - 1 ? std::numeric_limits<double>::infinity() : low[i] - high[i + 1];
    double min_gap = std::min(low_gap, high_gap);
    if (abs(min_gap / ((high[i] + low[i]) * 0.5)) > min_rel_sep) {
      twist_idx = get_twisted_factorization(d, l, eigenvals[i] - shift0, l_plus, u_minus);
      cout << "\t\t" << i << " UNSHIFTED gap: " << min_gap / (eigenvals[i] - shift0) << endl;
    }
    else if (abs(min_gap / ((high[i] + low[i]) * 0.5 - shift)) > min_rel_sep && min_element_growth < max_ele_growth) {
      double shifted_low = low[i] - shift;
      double shifted_high = high[i] - shift;
      eigenvalBisectRefine(d2, l2, shifted_low, shifted_high, i);
      double shifted_eigenval = (shifted_low + shifted_high) * 0.5;
      twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
      cout << "\t\t" << i << " prev shift gap: " << min_gap / ((high[i] + low[i]) * 0.5 - shift) << endl;
    }
    else {
      double max_shift = min_gap / min_rel_sep;
      bool need_aditional_step=false;
      if(max_shift<=0){
        max_shift=(high[i]-low[i])*10;
        need_aditional_step=true;
      }
      findShift(l, d, low[i], high[i], max_ele_growth, max_shift, l2, d2, shift, min_element_growth);

      //shift and refine eigenvalue
      double shifted_low = low[i] - shift;
      double shifted_high = high[i] - shift;
      eigenvalBisectRefine(d2, l2, shifted_low, shifted_high, i);
      if(need_aditional_step){
        cout << "aditional step!" << endl;
        low_gap = std::numeric_limits<double>::infinity();
        high_gap = std::numeric_limits<double>::infinity();
        if(i!=0){
          double shifted_prev_low = low[i-1] - shift;
          double shifted_prev_high = high[i-1] - shift;
          eigenvalBisectRefine(d2, l2, shifted_prev_low, shifted_prev_high, i-1);
          low_gap = shifted_prev_low - shifted_high;
        }
        if(i!=n-1){
          double shifted_next_low =low[i+1] - shift;
          double shifted_next_high = high[i+1] - shift;
          eigenvalBisectRefine(d2, l2, shifted_next_low, shifted_next_high, i+1);
          high_gap = shifted_low - shifted_next_high;
        }
        min_gap = std::min(low_gap, high_gap);
        if(min_gap<0){//TODO perturb
          cout << "################# TODO: need perturb";
        }
        max_shift = min_gap / min_rel_sep;
        double shift2;
        findShift(l2, d2, shifted_low, shifted_high, max_ele_growth, max_shift, l3, d3, shift2, min_element_growth);
        shifted_low -= shift2;
        shifted_high -= shift2;
        eigenvalBisectRefine(d3, l3, shifted_low, shifted_high, i);
        double shifted_eigenval = (shifted_low + shifted_high) * 0.5;
        cout << "\t\t" << i << " shifted gap: " << min_gap / shifted_eigenval << endl;
        twist_idx = get_twisted_factorization(d3, l3, shifted_eigenval, l_plus, u_minus);
      }
      else {
        double shifted_eigenval = (shifted_low + shifted_high) * 0.5;
        cout << "\t\t" << i << " shifted gap: " << min_gap / shifted_eigenval << endl;
        twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
      }
    }
    //calculate eigenvector
    calculateEigenvector(l_plus, u_minus, subdiag, i, twist_idx, eigenvecs);
  }
}
/*
Eigen::Array<int, BISECT_K, 1> getSturmCountLdlVec(const Eigen::VectorXd& d, const Eigen::VectorXd& l, Eigen::Array<double, BISECT_K, 1> shifts, int n_valid){
  int n = l.size();
  Eigen::Array<double, BISECT_K, 1> s;
  s.head(n_valid) = -shifts.head(n_valid);
  Eigen::Array<double, BISECT_K, 1> d_plus, l_plus;
  Eigen::Array<bool, BISECT_K, 1> cond;
  Eigen::Array<int, BISECT_K, 1> counts;
  counts.head(n_valid) = Eigen::ArrayXd(n_valid,0);
  for (int i = 0; i < n; i++) {
    d_plus.head(n_valid) = s.head(n_valid) + d[i];
    counts.head(n_valid) += (d_plus.head(n_valid) >= 0).cast<int>();
    cond.head(n_valid) = d_plus.head(n_valid).unaryExpr([](double x){return is_inf(x);}) &&
                              s.head(n_valid).unaryExpr([](double x){return is_inf(x);});
    s.head(n_valid) = l[i] * l[i] * s.head(n_valid) * (d[i] / d_plus.head(n_valid)) - shifts.head(n_valid);
    for(int j=0;j<n_valid;j++){
      if(cond[j]){
        s[j] = l[i] * l[i] * d[i] - shifts[j];
      }
    }
  }
  return Eigen::Array<int, BISECT_K, 1>;
}

void mrrr2(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs, double min_rel_sep = 1e-2, double max_ele_growth = 2) {
  int n = diag.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  getGresgorin(diag, subdiag, min_eigval, max_eigval);
  double shift0 = min_eigval - (max_eigval - min_eigval) * 0.1;
  Eigen::VectorXd l(n - 1), d(n), l0(n - 1), d0(n);
  get_ldl(diag, subdiag, shift0, l0, d0);
  for (int i = 0; i < n; i++) {
    if (i != n - 1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i] = d0[i] * get_random_perturbation_multiplier();
  }
  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();
  eigenvalsBisect4(diag, subdiagSquared, min_eigval, max_eigval);
}
*/

void symmetricEigenSolver(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs) {
  Eigen::MatrixXd packed;
  auto start = std::chrono::steady_clock::now();
  block_householder_tridiag3(A, packed);

  cout << "tridiag: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  Eigen::VectorXd diag = packed.diagonal();
  Eigen::VectorXd subdiag = packed.diagonal(1);
  cout << "extract: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  mrrr(diag, subdiag, eigenvals, eigenvecs);

  cout << "mrrr: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  start = std::chrono::steady_clock::now();
  block_apply_packed_Q3(packed, eigenvecs);
  cout << "apply q: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
}


}
}

#endif //STAN_OPENCL
#endif //STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP
