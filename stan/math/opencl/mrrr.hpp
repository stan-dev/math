#ifndef STAN_MATH_OPENCL_MRRR_HPP
#define STAN_MATH_OPENCL_MRRR_HPP

#ifdef STAN_OPENCL

#include <queue>

#include <stan/math/prim/mat/fun/mrrr.hpp>

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/subtract.hpp>
#include <stan/math/opencl/add.hpp>
#include <stan/math/opencl/transpose.hpp>

#include <stan/math/opencl/kernels/mrrr.hpp>

namespace stan {
namespace math {

/**
 * Calculates eigenvalues and eigenvectors of a (preferrably irreducible) tridiagonal matrix T using MRRR algorithm.
 * @param diag Diagonal of of T.
 * @param subdiag Subdiagonal of T.
 * @param eigenvals[out] Eigenvlues.
 * @param eigenvecs[out] Eigenvectors.
 * @param min_rel_sep Minimal relative separation of eigenvalues before computing eigenvectors.
 * @param max_ele_growth Maximal desired element growth of LDL decompositions.
 */
void mrrr_cl(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::VectorXd& subdiag, Eigen::Ref<Eigen::VectorXd> eigenvals,
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

  matrix_cl l_gpu(l);
  matrix_cl d_gpu(d);
  matrix_cl low_gpu(n, 1);
  matrix_cl high_gpu(n, 1);
  try {
    opencl_kernels::eigenvals_bisect(
            cl::NDRange(n),
            l_gpu.buffer(), d_gpu.buffer(), low_gpu.buffer(), high_gpu.buffer(),
            min_eigval - shift0, max_eigval - shift0, n);
  }
  catch (cl::Error& e) {
    check_opencl_error("block_apply_packed_Q_cl3", e);
  }
  copy(low, low_gpu);
  copy(high, high_gpu);
  eigenvals = (high + low).array() * 0.5 + shift0;
  std::queue<mrrrTask> blockQueue;
  blockQueue.push(mrrrTask{0, n, shift0, std::move(l), std::move(d), 0});
  l.resize(n - 1); //after move out
  d.resize(n);
  Eigen::MatrixXd l_big(n-1,n), d_big(n,n);
  Eigen::VectorXd min_gap_big(n);
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
        int twist_idx;
        double low_gap = i == block.start ? std::numeric_limits<double>::infinity() : low[i - 1] - high[i];
        double high_gap = i == block.end - 1 ? std::numeric_limits<double>::infinity() : low[i] - high[i + 1];
        double min_gap = std::min(low_gap, high_gap);
        min_gap_big[i]=min_gap;
        l_big.col(i) = block.l;
        d_big.col(i) = block.d;
      }
    }
  }
  matrix_cl subdiag_gpu(subdiag);
  matrix_cl l_big_gpu(l_big);
  matrix_cl l_big_t_gpu = transpose(l_big_gpu);
  matrix_cl d_big_gpu(d_big);
  matrix_cl d_big_t_gpu = transpose(d_big_gpu);
  copy(low_gpu,low);
  copy(high_gpu, high);
  matrix_cl min_gap_gpu(min_gap_big);

  matrix_cl temp1(n,n);
  matrix_cl temp2(n,n);
  matrix_cl temp3(n,n);
  matrix_cl temp4(n,n);
  matrix_cl temp5(n,n);
  matrix_cl eigenvecs_t_gpu(n,n);
  try {
    opencl_kernels::get_eigenvectors(
            cl::NDRange(n),
            subdiag_gpu.buffer(), l_big_t_gpu.buffer(), d_big_t_gpu.buffer(),
            low_gpu.buffer(), high_gpu.buffer(), min_gap_gpu.buffer(),
            temp1.buffer(),temp2.buffer(),temp3.buffer(),temp4.buffer(),temp5.buffer(),
            eigenvecs_t_gpu.buffer(), min_rel_sep, max_ele_growth);
  }
  catch (cl::Error& e) {
    check_opencl_error("block_apply_packed_Q_cl3", e);
  }
  matrix_cl eigenvecs_gpu = transpose(eigenvecs_t_gpu);
  Eigen::MatrixXd eigenvecs_tmp(n,n);
  copy(eigenvecs_tmp, eigenvecs_gpu);
  eigenvecs = eigenvecs_tmp;
}

/**
* Calculates eigenvalues and eigenvectors of a tridiagonal matrix T using MRRR algorithm on GPU.
* If a subdiagonal element is close to zero compared to neighbors on diagonal the problem can be split into smaller ones.
* @param diag Diagonal of of T.
* @param subdiag Subdiagonal of T.
* @param eigenvals[out] Eigenvlues.
* @param eigenvecs[out] Eigenvectors.
* @param splitThreshold Threshold for splitting the problem
*/
void tridiagonal_eigensolver_cl(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs, double splitThreshold = 1e-12) {
  int n = diag.size();
  eigenvecs.resize(n, n);
  eigenvals.resize(n);
  int last = 0;
  for (int i = 0; i < subdiag.size(); i++) {
    if (abs(subdiag[i] / diag[i]) < splitThreshold && abs(subdiag[i] / diag[i + 1]) < splitThreshold) {
      eigenvecs.block(last, i + 1, i + 1 - last, n - i - 1) = Eigen::MatrixXd::Constant(i + 1 - last, n - i - 1, 0);
      eigenvecs.block(i + 1, last, n - i - 1, i + 1 - last) = Eigen::MatrixXd::Constant(n - i - 1, i + 1 - last, 0);
      if (last == i) {
        eigenvecs(last, last) = 1;
        eigenvals[last] = diag[last];
      }
      else {
        mrrr_cl(diag.segment(last, i + 1 - last),
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
    mrrr_cl(diag.segment(last, n - last),
            subdiag.segment(last, subdiag.size() - last),
            eigenvals.segment(last, n - last),
            eigenvecs.block(last, last, n - last, n - last));
  }
}

}
}

#endif
#endif
