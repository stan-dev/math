#ifndef STAN_MATH_OPENCL_MRRR_HPP
#define STAN_MATH_OPENCL_MRRR_HPP

#ifdef STAN_OPENCL

#include <cmath>
#include <queue>

#include <stan/math/prim/mat/fun/mrrr.hpp>

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/subtract.hpp>
#include <stan/math/opencl/add.hpp>
#include <stan/math/opencl/transpose.hpp>
#include <stan/math/opencl/copy.hpp>

#include <stan/math/opencl/kernels/mrrr.hpp>

namespace stan {
namespace math {
namespace internal {

/**
 * Calculates eigenvalues and eigenvectors of a (preferrably irreducible) tridiagonal matrix T using MRRR algorithm.
 * @param diagonal Diagonal of of T.
 * @param subdiagonal Subdiagonal of T.
 * @param[out] eigenvalues Eigenvlues.
 * @param[out] eigenvectors Eigenvectors.
 * @param min_rel_sep Minimal relative separation of eigenvalues before computing eigenvectors.
 * @param max_ele_growth Maximal desired element growth of LDL decompositions.
 */
void mrrr_cl(const Eigen::Ref<const Eigen::VectorXd> diagonal, const Eigen::VectorXd& subdiagonal, Eigen::Ref<Eigen::VectorXd> eigenvalues,
             Eigen::Ref<Eigen::MatrixXd> eigenvectors, const double min_rel_sep = 1e-1, const double max_ele_growth = 2) {
  using std::copysign;
  const double shift_error = 1e-14;
  const int n = diagonal.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  get_gresgorin(diagonal, subdiagonal, min_eigval, max_eigval);
  Eigen::VectorXd l(n - 1), d(n), l0(n - 1), d0(n);
  const double shift0 = find_initial_shift(diagonal, subdiagonal, l0, d0, min_eigval, max_eigval, max_ele_growth);
  for (int i = 0; i < n; i++) {
    if (i != n - 1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i] = d0[i] * get_random_perturbation_multiplier();
  }

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
  eigenvalues = (high + low).array() * 0.5 + shift0;
  std::queue<mrrr_task> block_queue;
  block_queue.push(mrrr_task{0, n, shift0, std::move(l), std::move(d), 0});
  l.resize(n - 1); //after move out
  d.resize(n);
  Eigen::MatrixXd l_big(n - 1, n), d_big(n, n);
  Eigen::VectorXd min_gap_big(n);
  while (!block_queue.empty()) {
    mrrr_task block = block_queue.front();
    block_queue.pop();
    double shift = std::numeric_limits<double>::infinity();
    double min_element_growth = std::numeric_limits<double>::infinity();
    Eigen::VectorXd l2(n - 1), d2(n), l_plus(n - 1), u_minus(n - 1);
    for (int i = block.start; i < block.end; i++) {
      int cluster_end;
      for (cluster_end = i + 1; cluster_end < block.end; cluster_end++) {
        const int prev = cluster_end - 1;
        const double end_threshold = low[prev] * (1 - copysign(shift_error, low[prev]));
        if (high[cluster_end] < end_threshold) {
          break;
        }
      }
      cluster_end--; //now this is the index of the last element of the cluster
      if (cluster_end - i > 0) {//cluster
        const double max_shift = (high[i] - low[cluster_end]) * 10;
        double next_shift, min_ele_growth;
        find_shift(block.l, block.d, low[cluster_end], high[i], max_ele_growth, max_shift, l, d, next_shift, min_ele_growth);
        for (int j = i; j <= cluster_end; j++) {
          low[j] = low[j] * (1 - copysign(shift_error, low[j])) - next_shift;
          high[j] = high[j] * (1 + copysign(shift_error, high[j])) - next_shift;
          eigenval_bisect_refine(l, d, low[j], high[j], j);
        }
        block_queue.push(mrrr_task{i, cluster_end + 1, block.shift + next_shift, std::move(l), std::move(d), block.level + 1});
        l.resize(n - 1); //after move out
        d.resize(n);

        i = cluster_end;
      }
      else { //isolated eigenvalue
        int twist_idx;
        const double low_gap = i == block.start ? std::numeric_limits<double>::infinity() : low[i - 1] - high[i];
        const double high_gap = i == block.end - 1 ? std::numeric_limits<double>::infinity() : low[i] - high[i + 1];
        const double min_gap = std::min(low_gap, high_gap);
        min_gap_big[i] = min_gap;
        l_big.col(i) = block.l;
        d_big.col(i) = block.d;
      }
    }
  }
  matrix_cl subdiagonal_gpu(subdiagonal);
  matrix_cl l_big_gpu(l_big);
  matrix_cl l_big_t_gpu = transpose(l_big_gpu);
  matrix_cl d_big_gpu(d_big);
  matrix_cl d_big_t_gpu = transpose(d_big_gpu);
  copy(low_gpu, low);
  copy(high_gpu, high);
  matrix_cl min_gap_gpu(min_gap_big);

  matrix_cl temp1(n, n);
  matrix_cl temp2(n, n);
  matrix_cl temp3(n, n);
  matrix_cl temp4(n, n);
  matrix_cl temp5(n, n);
  matrix_cl eigenvectors_t_gpu(n, n);
  try {
    opencl_kernels::get_eigenvectors(
            cl::NDRange(n),
            subdiagonal_gpu.buffer(), l_big_t_gpu.buffer(), d_big_t_gpu.buffer(),
            low_gpu.buffer(), high_gpu.buffer(), min_gap_gpu.buffer(),
            temp1.buffer(), temp2.buffer(), temp3.buffer(), temp4.buffer(), temp5.buffer(),
            eigenvectors_t_gpu.buffer(), min_rel_sep, max_ele_growth);
  }
  catch (cl::Error& e) {
    check_opencl_error("block_apply_packed_Q_cl3", e);
  }
  matrix_cl eigenvectors_gpu = transpose(eigenvectors_t_gpu);
  Eigen::MatrixXd eigenvectors_tmp(n, n);
  copy(eigenvectors_tmp, eigenvectors_gpu);
  eigenvectors = eigenvectors_tmp;
}

/**
* Calculates eigenvalues and eigenvectors of a tridiagonal matrix T using MRRR algorithm on GPU.
* If a subdiagonal element is close to zero compared to neighbors on diagonal the problem can be split into smaller ones.
* @param diagonal Diagonal of of T.
* @param subdiagonal Subdiagonal of T.
* @param[out] eigenvalues Eigenvlues.
* @param[out] eigenvectors Eigenvectors.
* @param split_threshold Threshold for splitting the problem
*/
void tridiagonal_eigensolver_cl(const Eigen::VectorXd& diagonal, const Eigen::VectorXd& subdiagonal, Eigen::VectorXd& eigenvalues, Eigen::MatrixXd& eigenvectors, const double split_threshold = 1e-12) {
  const int n = diagonal.size();
  eigenvectors.resize(n, n);
  eigenvalues.resize(n);
  int last = 0;
  for (int i = 0; i < subdiagonal.size(); i++) {
    if (std::fabs(subdiagonal[i] / diagonal[i]) < split_threshold && std::fabs(subdiagonal[i] / diagonal[i + 1]) < split_threshold) {
      eigenvectors.block(last, i + 1, i + 1 - last, n - i - 1) = Eigen::MatrixXd::Constant(i + 1 - last, n - i - 1, 0);
      eigenvectors.block(i + 1, last, n - i - 1, i + 1 - last) = Eigen::MatrixXd::Constant(n - i - 1, i + 1 - last, 0);
      if (last == i) {
        eigenvectors(last, last) = 1;
        eigenvalues[last] = diagonal[last];
      }
      else {
        mrrr_cl(diagonal.segment(last, i + 1 - last),
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
  }
  else {
    mrrr_cl(diagonal.segment(last, n - last),
            subdiagonal.segment(last, subdiagonal.size() - last),
            eigenvalues.segment(last, n - last),
            eigenvectors.block(last, last, n - last, n - last));
  }
}

}
}
}

#endif
#endif
