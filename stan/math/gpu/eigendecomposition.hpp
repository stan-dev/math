//
// Created by tadej on 7. 12. 2018.
//

#include <Eigen/QR>
#include <iostream>

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
  std::cout << "(" << a.rows() << ", " << a.cols() << ")"  << std::endl;
}

void s(const Eigen::RowVectorXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << ")"  << std::endl;
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
      if(j!=0) {
        auto householder_whole = packed.col(k + j).tail(packed.rows() - k - j);
        householder_whole -= U.block(j - 1,0,householder_whole.size(),j) * V.block(j-1,0,1,j).transpose() +
                             V.block(j - 1,0,householder_whole.size(),j) * U.block(j-1,0,1,j).transpose();
      }
      double q = householder.squaredNorm();
      double alpha = -copysign(sqrt(q), packed(k+j,k+j));
      q -= householder[0] * householder[0];
      householder[0] -= alpha;
      q+=householder[0]*householder[0];
      q=sqrt(q);
      householder *= SQRT_2 / q;

      auto& u = householder;
      Eigen::VectorXd v(householder.size()+1);
      v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - j - 1, packed.cols() - k - j -1).selfadjointView<Eigen::Lower>() * u
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
    packed.block(k + actual_r , k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -=
            (Y_partial_update + Y_partial_update.transpose()).bottomRightCorner(Y_partial_update.rows()-actual_r + 1, Y_partial_update.cols()-actual_r + 1);
  }
  packed(packed.rows()-2, packed.cols()-1) = packed(packed.rows()-1, packed.cols()-2);
}

void block_apply_packed_Q3(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100){
  //if input A==Identity, constructs Q
  Eigen::MatrixXd scratchSpace(A.rows(), r);
  for (int k = (packed.rows() - 3)/r*r; k >=0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k -1, actual_r);
    W.col(0) = packed.col(k).tail(W.rows());
    for (size_t j = 1; j < actual_r; j++) {
      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
      W.col(j).noalias() = - W.leftCols(j) * scratchSpace.col(0).head(j);
      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
    }
    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
  }
}

const double perturbation_range = 1e-13;
inline double get_random_perturbation_multiplier(){
  static const double rand_norm = perturbation_range / RAND_MAX;
  static const double almost_one = 1 - perturbation_range * 0.5;
  return almost_one + std::rand() * rand_norm;
}

void get_ldl(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, double shift, Eigen::VectorXd& l, Eigen::VectorXd& d_plus){
  d_plus[0] = diag[0] - shift;
  for(int i = 0; i < subdiag.size(); i++){
    l[i] = subdiag[i] / d_plus[i];
    //d_plus[i] *= get_random_perturbation_multiplier();
    d_plus[i+1] = diag[i+1] - shift - l[i] * subdiag[i];
    //l[i] *= get_random_perturbation_multiplier();
  }
  //d_plus[subdiag.size()] *= get_random_perturbation_multiplier();
}

//stationary qds
double get_perturbed_shifted_ldl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& d_plus){
  int n = l.size();
  double s = -shift;
  double element_growth = 0;
  double element_growth_denominator = 0;
  for(int i=0; i < n; i++){
    double di_perturbed = d[i];// * get_random_perturbation_multiplier();
    double li_perturbed = l[i];// * get_random_perturbation_multiplier();
    d_plus[i] = s + di_perturbed;
    element_growth += abs(d_plus[i]);
    element_growth_denominator += d_plus[i];
    l_plus[i] = li_perturbed * (di_perturbed / d_plus[i]);
    if(is_inf(d_plus[i]) && is_inf(s)) { // this happens if d_plus[i]==0 -> in next iteration d_plus==inf and s==inf
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
int get_twisted_factorization(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& u_minus){
  int n = l.size();
  //calculate shifted ldl
  Eigen::VectorXd s(n+1);
  s[0] = -shift;
  for(int i=0; i < n; i++){
    double d_plus = s[i] + d[i];
    l_plus[i] = l[i] * (d[i] / d_plus);
    if(is_nan(l_plus[i])){ //d_plus==0
      //one (or both) of d[i], l[i] is very close to 0
      if(abs(l[i])<abs(d[i])){
        l_plus[i] = d[i] * copysign(1.,l[i]) * copysign(1.,d_plus);
      }
      else{
        l_plus[i] = l[i] * copysign(1.,d[i]) * copysign(1.,d_plus);
      }
    }
    s[i+1] = l_plus[i] * l[i] * s[i] - shift;
    if(is_nan(s[i+1])){
      if(abs(l_plus[i])>abs(s[i])){ //l_plus[i]==inf
        if(abs(s[i])>abs(l[i])){ //l[i]==0
          s[i+1] = s[i] * copysign(1.,l[i]) * copysign(1.,l_plus[i]) - shift;
        }
        else{ //s[i]==0
          s[i+1] = l[i] * copysign(1.,s[i]) * copysign(1.,l_plus[i]) - shift;
        }
      }
      else{ //s[i]==inf
        if(abs(l_plus[i])>abs(l[i])){ //l[i]==0
          s[i+1] = l_plus[i] * copysign(1.,l[i]) * copysign(1.,s[i]) - shift;
        }
        else{ //l_plus[i]==0
          s[i+1] = l[i] * copysign(1.,s[i]) * copysign(1.,l_plus[i]) - shift;
        }
      }
    }
  }
  //calculate shifted udu and twist index
  double p=d[n]-shift;
  double min_gamma = abs(s[n] + d[n]);
  int twist_index=n;

  for(int i = n - 1; i >= 0; i--){
    double d_minus = d[i] * l[i] * l[i] + p;
    double t = d[i] / d_minus;
    u_minus[i] = l[i] * t;
    if(is_nan(u_minus[i])){
      if(is_nan(t)){
        t = copysign(1.,d[i]) * copysign(1.,d_minus);
        u_minus[i] = l[i] * t;
      }
      else{ //t==inf, l[i]==0
        u_minus[i] = d[i] * copysign(1.,l[i]) * copysign(1.,t);
      }
    }
    double gamma = abs(s[i] + t * p);//TODO: all inf/nan -> need other shift!
    if(is_nan(gamma)){ //t==inf, p==0 OR t==0, p==inf
      double d_sign = d[i] * copysign(1.,d_minus) * copysign(1.,t);
      gamma = abs(s[i] + d_sign);
      p = d_sign - shift;
    }
    else{ //usual case
      p = p * t - shift;
    }
    if(gamma < min_gamma){
      min_gamma = gamma;
      twist_index = i;
    }
  }
  return twist_index;
}

int getSturmCountLdlLdl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift){
  int n = l.size();
  double s = -shift;
  double l_plus, d_plus;
  int count = 0;
  for(int i=0; i < n; i++){
    d_plus = s + d[i];
    count += d_plus>=0;
    if (is_inf(d_plus) && is_inf(s)){ // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
      s = l[i] * l[i] * d[i] - shift;
    }
    else{
      s = l[i] * l[i] * s * (d[i] / d_plus) - shift;
    }
  }
  d_plus = s + d[n];
  count += d_plus>=0;
  return count;
}

double eigenvalBisectRefine(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double& low, double& high, int i){
  int n = d.size();
  double eps = 3e-16;
  while(abs((high-low)/(high+low))>eps && abs(high-low)>std::numeric_limits<double>::min()){ // second term is for the case where the eigenvalue is 0 and division yields NaN
    double mid=(high+low)*0.5;
    if(getSturmCountLdlLdl(d, l, mid)>i){
      low=mid;
    }
    else{
      high=mid;
    }
  }
}

void getGresgorin(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, double& min_eigval, double& max_eigval) {
  int n = diag.size();
  min_eigval = diag[0] - abs(subdiag[0]);
  max_eigval = diag[0] + abs(subdiag[0]);
  for(int i=1;i<n-1;i++){
    min_eigval = std::min(min_eigval, diag[i] - abs(subdiag[i]) - abs(subdiag[i - 1]));
    max_eigval = std::max(max_eigval, diag[i] + abs(subdiag[i]) + abs(subdiag[i - 1]));
  }
  min_eigval = std::min(min_eigval, diag[n - 1] - abs(subdiag[n - 2]));
  max_eigval = std::max(max_eigval, diag[n - 1] + abs(subdiag[n - 2]));
}

#define BISECT_K 8
void eigenvalsBisect3(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiagSquared, double min_eigval, double max_eigval, Eigen::VectorXd& low, Eigen::VectorXd& high){
  int n = diag.size();
  double eps = 3e-16;

  for(int i=0; i<n; i += BISECT_K){
    int n_valid = std::min(BISECT_K,n-i);
    auto low_work = low.segment(i,n_valid).array();
    auto high_work = high.segment(i,n_valid).array();
    low_work=Eigen::Array<double, Eigen::Dynamic, 1>::Constant(n_valid, min_eigval);
    high_work=Eigen::Array<double, Eigen::Dynamic, 1>::Constant(n_valid, max_eigval);
    while((abs((high_work-low_work)/low_work)>eps && abs(high_work-low_work)>std::numeric_limits<double>::min()).any()){
      Eigen::Array<double, BISECT_K, 1> shifts;
      shifts.head(n_valid) = (high_work+low_work)*0.5;
      Eigen::Array<double, BISECT_K, 1> d;
      d.head(n_valid) = diag[0] - shifts.head(n_valid);
      Eigen::Array<int, BISECT_K, 1> counts;
      counts.head(n_valid) = (d.head(n_valid)>=0).cast<int>();
      for(int j = 1; j < diag.size(); j++){ //Sturm count via LDL factorization
        d.head(n_valid) = diag[j] - shifts.head(n_valid) - subdiagSquared[j-1] / d.head(n_valid);
        counts.head(n_valid) += (d.head(n_valid)>=0).cast<int>();
      }
      for(int k=0;k<n_valid;k++){
        if(counts[k]>i+k){
          low_work[k]=shifts[k];
        }
        else{
          high_work[k]=shifts[k];
        }
      }
    }
  }
}

#undef BISECT_K


void mrrr(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag,  Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs){
  double min_rel_sep=1e-2;
  double max_ele_growth = 2;
  int n=diag.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  getGresgorin(diag, subdiag, min_eigval, max_eigval);
  double span = max_eigval - min_eigval;
  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();
  eigenvalsBisect3(diag, subdiagSquared, min_eigval, max_eigval, low, high);
  Eigen::VectorXd l(n-1), d(n), l0(n-1), d0(n);
  double shift0 = min_eigval - (max_eigval - min_eigval) * 0.1;
  get_ldl(diag, subdiag, shift0, l0, d0);
  Eigen::VectorXd l2(n-1), d2(n), l3(n-1), d3(n), l_plus(n-1), d_plus(n), u_minus(n-1), d_minus(n), s(n), p(n);

  double T_norm = sqrt((diag.array()*diag.array()).sum() + subdiagSquared.sum());
  for(int i=0;i<n;i++) {
    eigenvals[i]=(high[i]+low[i])*0.5;
    //from here on eigenvals[i] constins original eigenval
    if(i!=n-1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i]=d0[i]*get_random_perturbation_multiplier();
    low[i] -= shift0;
    high[i] -= shift0;
    //from here on high[i] and low[i] contain bounds on shifted representation
  }
  low[0] = low[0] * (1 - copysign(perturbation_range * n, low[0]));
  high[0] = high[0] * (1 + copysign(perturbation_range * n, high[0]));
  eigenvalBisectRefine(d, l, low[0], high[0], 0);
  
  double shift = std::numeric_limits<double>::infinity();
  double min_element_growth = std::numeric_limits<double>::infinity();
  int cluster_start=-1, cluster_end=-1;
  for(int i=0;i<n;i++){
    if(i>cluster_end){
      //double start_threshold = high[i] * (1 + copysign(sqrt(perturbation_range) * n, high[i]));
      double end_threshold = low[i] * (1 - copysign(sqrt(perturbation_range) * n, low[i]));
      /*for(cluster_start=i;cluster_start>1;cluster_start--){
        if(low[cluster_start-1]>start_threshold){
          break;
        }
      }*/
      cluster_start=i;
      for(cluster_end=i;cluster_end<n-1;cluster_end++){
        int next = cluster_end + 1;
        low[next] = low[next] * (1 - copysign(perturbation_range * n, low[next]));
        high[next] = high[next] * (1 + copysign(perturbation_range * n, high[next]));
        eigenvalBisectRefine(d, l, low[next], high[next], next);
        if(high[next] < end_threshold){
          break;
        }
      }
      while(true) { // perturb the representation if any eigenvalues in the cluster are too close
        bool is_ok = true;
        for(int j=cluster_start;j<cluster_end;j++){
          if(high[j+1] >= low[j]){//problem - need new perturb
            is_ok = false;
            break;
          }
        }
        if(is_ok){
          break;
        }
        cout << "perturbing at i = " << i << " cluster: " << cluster_start << ", " << cluster_end << endl;
        for (int j = 0; j < n - 1; j++) {
          l[j] = l0[j] * get_random_perturbation_multiplier();
          d[j] = d0[j] * get_random_perturbation_multiplier();
        }
        d[n - 1] = d0[n - 1] * get_random_perturbation_multiplier();
        //eigenvalBisectRefine(d, l, low[i], high[i], i);
        for(int j=cluster_start;j<=cluster_end;j++){
          low[j] = low[j] * (1 - copysign(perturbation_range * n, low[j]));
          high[j] = high[j] * (1 + copysign(perturbation_range * n, high[j]));
          eigenvalBisectRefine(d, l, low[j], high[j], j);
        }
        /*if(i!=0 && high[i]>=low[i-1]){//problem - need new perturb
          continue;
        }*/
      }
      if(cluster_end!=n-1){
        int next = cluster_end + 1;
        low[next] = low[next] * (1 - copysign(perturbation_range * n, low[next]));
        high[next] = high[next] * (1 + copysign(perturbation_range * n, high[next]));
        eigenvalBisectRefine(d, l, low[next], high[next], next);
      }
//      low_gap = i==0 ? std::numeric_limits<double>::infinity() : low[i-1] - high[i];
//      high_gap = i==n-1 ? std::numeric_limits<double>::infinity() : low[i] - high[i+1];
//      min_gap = std::min(low_gap, high_gap);
      min_element_growth = std::numeric_limits<double>::infinity();
    }
    int twist_idx;
    double low_gap = i==0 ? std::numeric_limits<double>::infinity() : low[i-1] - high[i];
    double high_gap = i==n-1 ? std::numeric_limits<double>::infinity() : low[i] - high[i+1];
    double min_gap = std::min(low_gap, high_gap);
    //if(min_gap<0){ //this and the next eigenvalue are the same!

    //}
    if(abs(min_gap / ((high[i]+low[i])*0.5)) > min_rel_sep){
      twist_idx = get_twisted_factorization(d, l, eigenvals[i] - shift0, l_plus, u_minus);
      cout << "\t\t" << i << " UNSHIFTED gap: " << min_gap/(eigenvals[i]-shift0) << endl;
    }
    else if(abs(min_gap / ((high[i]+low[i])*0.5 - shift)) > min_rel_sep && min_element_growth < max_ele_growth){
      double shifted_low=low[i] - shift;
      double shifted_high=high[i] - shift;
      eigenvalBisectRefine(d2, l2, shifted_low, shifted_high, i);
      double shifted_eigenval = (shifted_low + shifted_high) * 0.5;
      twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
      cout << "\t\t" << i << " prev shift gap: " << min_gap/((high[i]+low[i])*0.5-shift) << endl;
    }
    else {
      //TODO: do we need more shift options?
      std::vector<double> shifts;
      shifts.push_back(low[i]);
      shifts.push_back(high[i]);
      double max_shift = min_gap / min_rel_sep;
      if(i!=0){
        if((high[i-1] - low[i]) * 0.5 < max_shift){
          shifts.push_back(high[i] + (high[i-1] - low[i]) * 0.5);
        }
        if((high[i-1] - low[i]) * 0.25 < max_shift){
          shifts.push_back(high[i] + (high[i-1] - low[i]) * 0.25);
        }
        if((high[i-1] - low[i]) < max_shift){
          shifts.push_back(high[i] + (high[i-1] - low[i]));
        }
        if((high[i-1] - low[i]) * 3 < max_shift){
          shifts.push_back(high[i] + (high[i-1] - low[i]) * 3);
        }
      }
      if(i!=n-1){
        if((high[i]-low[i+1]) * 0.5 < max_shift){
          shifts.push_back(low[i] - (high[i]-low[i+1]) * 0.5);
        }
        if((high[i]-low[i+1]) * 0.25 < max_shift){
          shifts.push_back(low[i] - (high[i]-low[i+1]) * 0.25);
        }
        if((high[i]-low[i+1]) < max_shift){
          shifts.push_back(low[i] - (high[i]-low[i+1]));
        }
        if((high[i]-low[i+1]) * 3 < max_shift){
          shifts.push_back(low[i] - (high[i]-low[i+1]) * 3);
        }
      }
      shifts.push_back(high[i] - max_shift);
      shifts.push_back(low[i] + max_shift);
      shifts.push_back(high[i] - max_shift * 0.5);
      shifts.push_back(low[i] + max_shift * 0.5);
      shifts.push_back(high[i] - max_shift * 0.25);
      shifts.push_back(low[i] + max_shift * 0.25);
      shifts.push_back(high[i] - max_shift * 0.75);
      shifts.push_back(low[i] + max_shift * 0.75);
      min_element_growth = std::numeric_limits<double>::infinity();
      for(double sh : shifts){
        //sh -= shift0;
        double element_growth = get_perturbed_shifted_ldl(d, l, sh, l3, d3);
        cout << i << " element growth: " << element_growth << " at " << sh + shift0 << endl;
        if(element_growth<min_element_growth){
          l2.swap(l3);
          d2.swap(d3);
          shift = sh;// + shift0;
          min_element_growth=element_growth;
          if(element_growth<=max_ele_growth){
            break;
          }
        }
      }
      cout << "\t\t" << i << " element growth: " << min_element_growth << " at " << shift << endl;

      //shift and refine eigenvalue
      double shifted_low=low[i] - shift;
      double shifted_high=high[i] - shift;
      eigenvalBisectRefine(d2, l2, shifted_low, shifted_high, i);
      double shifted_eigenval = (shifted_low + shifted_high) * 0.5;
      cout << "\t\t" << i << " shifted gap: " << min_gap/shifted_eigenval << endl;
      twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
    }
    //calculate eigenvector
    auto vec = eigenvecs.col(i);
    vec[twist_idx] = 1;
    for(int j=twist_idx+1;j<n;j++){
      if(vec[j-1]!=0) {
        vec[j] = -u_minus[j - 1] * vec[j - 1];
      }
      else{
        vec[j] = -subdiag[j - 2] * vec[j - 2] / subdiag[j - 1];
        if(is_nan(vec[j]) || is_inf(vec[j])){ //subdiag[j - 1]==0
          vec[j]=0;
        }
      }
    }
    for(int j = twist_idx - 1; j >= 0;j--){
      if(vec[j+1]!=0){
        vec[j] = -l_plus[j] * vec[j + 1];
      }
      else{
        vec[j] = -subdiag[j + 1] * vec[j + 2] / subdiag[j];
        if(is_nan(vec[j]) || is_inf(vec[j])) { //subdiag[j]==0
          vec[j]=0;
        }
      }
    }
    //vec /= vec.norm();
  }
}

void symmetricEigenSolver(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs){
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
