//
// Created by tadej on 7. 12. 2018.
//

#include <Eigen/QR>
#include <iostream>
#include <queue>

#include <stan/math/gpu/matrix_gpu.hpp>
#include <stan/math/gpu/multiply.hpp>
#include <stan/math/gpu/subtract.hpp>
#include <stan/math/gpu/add.hpp>
#include <stan/math/gpu/transpose.hpp>

#include <stan/math/gpu/kernels/eigendecomposition.hpp>

using namespace std;

#ifndef STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP
#define STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP

#ifdef STAN_OPENCL

#define TIME_IT
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

/**
 * Tridiagonalize a symmetric matrix using block Housholder algorithm. A = Q * T * Q^T, where T is tridiagonal and Q is orthonormal.
 * @param A Input matrix
 * @param[out] packed Packed form of the tridiagonal matrix. Elements of the resulting symmetric tridiagonal matrix T are in the diagonal and first superdiagonal.
 * Columns bellow diagonal contain householder vectors that can be used to construct orthogonal matrix Q.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
void block_householder_tridiag4(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
  packed = A;
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < packed.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd V(packed.rows() - k - 1, actual_r);
    Eigen::MatrixXd U(packed.rows() - k - 1, actual_r);
    V.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);
    U.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(U.rows(), U.cols(), 0);

    for (size_t j = 0; j < actual_r; j++) {
#ifdef TIME_IT
      start = std::chrono::steady_clock::now();
#endif
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
#ifdef TIME_IT
      t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      start = std::chrono::steady_clock::now();
#endif

      auto& u = householder;
      Eigen::VectorXd v(householder.size() + 1);
      v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - j - 1, packed.cols() - k - j - 1).selfadjointView<Eigen::Lower>() * u
                                   - U.bottomLeftCorner(u.size(), j) * (V.bottomLeftCorner(u.size(), j).transpose() * u)
                                   - V.bottomLeftCorner(u.size(), j) * (U.bottomLeftCorner(u.size(), j).transpose() * u);
      v[0] = q / SQRT_2;// - alpha * householder[1];
      double cnst = v.tail(householder.size()).transpose() * u;
      v.tail(householder.size()) -= 0.5 * cnst * u;
#ifdef TIME_IT
      t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      start = std::chrono::steady_clock::now();
#endif

      packed(k + j, k + j + 1) = packed(k + j + 1, k + j) * q / SQRT_2 + alpha - v[0] * u[0];
      U.col(j).tail(U.rows() - j) = u;
      V.col(j).tail(V.rows() - j) = v.tail(V.rows() - j);
#ifdef TIME_IT
      t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    Eigen::MatrixXd partial_update = U.bottomRows(U.rows() - actual_r + 1) * V.bottomRows(V.rows() - actual_r + 1).transpose();
    packed.block(k + actual_r, k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -= partial_update + partial_update.transpose();
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
  packed(packed.rows() - 2, packed.cols() - 1) = packed(packed.rows() - 1, packed.cols() - 2);
#ifdef TIME_IT
  std::cout << "block_householder_tridiag4 detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

/**
 * Tridiagonalize a symmetric matrix using block Housholder algorithm. A = Q * T * Q^T, where T is tridiagonal and Q is orthonormal.
 * @param A Input matrix
 * @param[out] packed Packed form of the tridiagonal matrix. Elements of the resulting symmetric tridiagonal matrix T are in the diagonal and first superdiagonal.
 * Columns bellow diagonal contain householder vectors that can be used to construct orthogonal matrix Q.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
void block_householder_tridiag5(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
  packed = A;
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < packed.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd V(packed.rows() - k - 1, actual_r);
    Eigen::MatrixXd U(packed.rows() - k - 1, actual_r);
    V.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);

    for (size_t j = 0; j < actual_r; j++) {
#ifdef TIME_IT
      start = std::chrono::steady_clock::now();
#endif
      auto householder = packed.col(k + j).tail(packed.rows() - k - j - 1);
      if (j != 0) {
        auto householder_whole = packed.col(k + j).tail(packed.rows() - k - j);
//        householder_whole -= U.block(j - 1, 0, householder_whole.size(), j) * V.block(j - 1, 0, 1, j).transpose() +
//                             V.block(j - 1, 0, householder_whole.size(), j) * U.block(j - 1, 0, 1, j).transpose();
        householder_whole -= packed.block(j + k, k, householder_whole.size(), j) * V.block(j - 1, 0, 1, j).transpose() +
                             V.block(j - 1, 0, householder_whole.size(), j) * packed.block(j + k, k, 1, j).transpose();
      }
      double q = householder.squaredNorm();
      double alpha = -copysign(sqrt(q), packed(k + j, k + j));
      q -= householder[0] * householder[0];
      householder[0] -= alpha;
      q += householder[0] * householder[0];
      q = sqrt(q);
      householder *= SQRT_2 / q;
#ifdef TIME_IT
      t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      start = std::chrono::steady_clock::now();
#endif

      auto& u = householder;
      Eigen::VectorXd v(householder.size() + 1);
//      Eigen::VectorXd R1 = V.bottomLeftCorner(u.size(), j).transpose() * u;
//      Eigen::VectorXd R2 = packed.block(k + j + 1, k, u.size(), j).transpose() * u;
      v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - j - 1, packed.cols() - k - j - 1).selfadjointView<Eigen::Lower>() * u
                                   - packed.block(k + j + 1, k, u.size(), j) * (V.bottomLeftCorner(u.size(), j).transpose() * u)
                                   - V.bottomLeftCorner(u.size(), j) * (packed.block(k + j + 1, k, u.size(), j).transpose() * u);
      v[0] = q / SQRT_2;// - alpha * householder[1];
      double cnst = v.tail(householder.size()).transpose() * u;
      v.tail(householder.size()) -= 0.5 * cnst * u;
#ifdef TIME_IT
      t2 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      start = std::chrono::steady_clock::now();
#endif

      packed(k + j, k + j + 1) = packed(k + j + 1, k + j) * q / SQRT_2 + alpha - v[0] * u[0];
      V.col(j).tail(V.rows() - j) = v.tail(V.rows() - j);
#ifdef TIME_IT
      t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    Eigen::MatrixXd partial_update = packed.block(k + actual_r, k,packed.rows() - k - actual_r, actual_r) * V.bottomRows(V.rows() - actual_r + 1).transpose();
    packed.block(k + actual_r, k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -= partial_update + partial_update.transpose();
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
  packed(packed.rows() - 2, packed.cols() - 1) = packed(packed.rows() - 1, packed.cols() - 2);
#ifdef TIME_IT
  std::cout << "block_householder_tridiag5 detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

/**
 * Tridiagonalize a symmetric matrix using block Housholder algorithm. A = Q * T * Q^T, where T is tridiagonal and Q is orthonormal.
 * @param A Input matrix
 * @param[out] packed Packed form of the tridiagonal matrix. Elements of the resulting symmetric tridiagonal matrix T are in the diagonal and first superdiagonal.
 * Columns bellow diagonal contain householder vectors that can be used to construct orthogonal matrix Q.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
void block_householder_tridiag_gpu(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
  packed = A;
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < packed.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd V(packed.rows() - k - 1, actual_r);
    Eigen::MatrixXd U(packed.rows() - k - 1, actual_r);
    V.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);
    U.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(U.rows(), U.cols(), 0);

    for (size_t j = 0; j < actual_r; j++) {
#ifdef TIME_IT
      start = std::chrono::steady_clock::now();
#endif
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
#ifdef TIME_IT
      t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      start = std::chrono::steady_clock::now();
#endif

      auto& u = householder;
      Eigen::VectorXd v(householder.size() + 1);
      v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - j - 1, packed.cols() - k - j - 1).selfadjointView<Eigen::Lower>() * u
                                   - U.bottomLeftCorner(u.size(), j) * (V.bottomLeftCorner(u.size(), j).transpose() * u)
                                   - V.bottomLeftCorner(u.size(), j) * (U.bottomLeftCorner(u.size(), j).transpose() * u);
      v[0] = q / SQRT_2;// - alpha * householder[1];
      double cnst = v.tail(householder.size()).transpose() * u;
      v.tail(householder.size()) -= 0.5 * cnst * u;
#ifdef TIME_IT
      t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      start = std::chrono::steady_clock::now();
#endif

      packed(k + j, k + j + 1) = packed(k + j + 1, k + j) * q / SQRT_2 + alpha - v[0] * u[0];
      U.col(j).tail(U.rows() - j) = u;
      V.col(j).tail(V.rows() - j) = v.tail(V.rows() - j);
#ifdef TIME_IT
      t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    matrix_gpu U_gpu(U.bottomRows(U.rows() - actual_r + 1).eval());
    matrix_gpu V_T_gpu(V.bottomRows(V.rows() - actual_r + 1).transpose().eval());
    matrix_gpu partial_update_gpu = U_gpu * V_T_gpu;
    Eigen::MatrixXd partial_update(partial_update_gpu.rows(), partial_update_gpu.cols());
    copy(partial_update, partial_update_gpu);
    packed.block(k + actual_r, k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -= partial_update + partial_update.transpose();
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
  packed(packed.rows() - 2, packed.cols() - 1) = packed(packed.rows() - 1, packed.cols() - 2);
#ifdef TIME_IT
  std::cout << "block_householder_tridiag_gpu detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

/**
 * Tridiagonalize a symmetric matrix using block Housholder algorithm. A = Q * T * Q^T, where T is tridiagonal and Q is orthonormal.
 * @param A Input matrix
 * @param[out] packed Packed form of the tridiagonal matrix. Elements of the resulting symmetric tridiagonal matrix T are in the diagonal and first superdiagonal.
 * Columns bellow diagonal contain householder vectors that can be used to construct orthogonal matrix Q.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
void block_householder_tridiag_gpu2(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
  matrix_gpu packed_gpu(A);
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0, t8=0, t9=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < A.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(A.rows() - k - 2)});

#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    matrix_gpu V_gpu(A.rows() - k - 1, actual_r);
    V_gpu.zeros<TriangularViewGPU::Upper>();
#ifdef TIME_IT
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif

    for (size_t j = 0; j < actual_r; j++) {
#ifdef TIME_IT
      start = std::chrono::steady_clock::now();
#endif
      matrix_gpu Uu(j,1), Vu(j,1), q_gpu(1,1);
      try{
//        opencl_kernels::eigendecomp_householder_v1(
//                cl::NDRange(128), cl::NDRange(128),
//                packed_gpu.buffer(), V_gpu.buffer(), q_gpu.buffer(), Uu.buffer(), Vu.buffer(),
//                packed_gpu.rows(), V_gpu.rows(), j, k);
        opencl_kernels::eigendecomp_householder(
                cl::NDRange(128), cl::NDRange(128),
                packed_gpu.buffer(), V_gpu.buffer(), q_gpu.buffer(),
                packed_gpu.rows(), V_gpu.rows(), j, k);
#ifdef TIME_IT
        t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
        if(j!=0) {
          opencl_kernels::eigendecomp_v1(
                  cl::NDRange(64 * j), cl::NDRange(64),
                  packed_gpu.buffer(), V_gpu.buffer(), Uu.buffer(), Vu.buffer(),
                  packed_gpu.rows(), V_gpu.rows(), k);
        }
#ifdef TIME_IT
        t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
        opencl_kernels::eigendecomp_v2(
                cl::NDRange((A.rows() - k - j - 1)*64),cl::NDRange(64),
                packed_gpu.buffer(), V_gpu.buffer(), Uu.buffer(), Vu.buffer(),
                packed_gpu.rows(), V_gpu.rows(), k, j);
//        opencl_kernels::eigendecomp_v2_2(
//                cl::NDRange((A.rows() - k - j - 1 + 31)/32*32),cl::NDRange(32),
//                packed_gpu.buffer(), V_gpu.buffer(), Uu.buffer(), Vu.buffer(),
//                packed_gpu.rows(), V_gpu.rows(), k, j);
#ifdef TIME_IT
        t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
        opencl_kernels::eigendecomp_v3(
                cl::NDRange(128),cl::NDRange(128),
                packed_gpu.buffer(), V_gpu.buffer(), q_gpu.buffer(),
                packed_gpu.rows(), V_gpu.rows(), k, j);
      }
      catch (cl::Error& e) {
        check_opencl_error("block_apply_packed_Q_gpu3", e);
      }
#ifdef TIME_IT
      t5+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    matrix_gpu U_gpu(V_gpu.rows() - actual_r + 1, V_gpu.cols());
    U_gpu.sub_block(packed_gpu, k + actual_r, k, 0, 0, A.rows() - k - actual_r, actual_r);
//    matrix_gpu V_T_gpu(V.bottomRows(V.rows() - actual_r + 1).transpose().eval());
    matrix_gpu Vb_gpu(V_gpu.rows() - actual_r + 1, V_gpu.cols());
    Vb_gpu.sub_block(V_gpu, actual_r - 1, 0, 0, 0, V_gpu.rows() - actual_r + 1, actual_r);
#ifdef TIME_IT
    t6+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
    matrix_gpu partial_update_gpu = U_gpu * transpose(Vb_gpu);
#ifdef TIME_IT
    t7+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
    try{
      opencl_kernels::subtract_twice(
              cl::NDRange(partial_update_gpu.rows(), partial_update_gpu.cols()),
              packed_gpu.buffer(), partial_update_gpu.buffer(),
              packed_gpu.rows(), partial_update_gpu.rows(), k + actual_r);
    }
    catch (cl::Error& e) {
      check_opencl_error("block_apply_packed_Q_gpu3", e);
    }
#ifdef TIME_IT
    t8+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
  }
#ifdef TIME_IT
  start = std::chrono::steady_clock::now();
#endif
  packed.resize(A.rows(),A.cols());
  copy(packed, packed_gpu);
  packed(packed.rows() - 2, packed.cols() - 1) = packed(packed.rows() - 1, packed.cols() - 2);
#ifdef TIME_IT
  t9+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  start = std::chrono::steady_clock::now();
#endif
#ifdef TIME_IT
  std::cout << "block_householder_tridiag_gpu2 detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << " " << t5/1000  << " " << t6/1000  << " " << t7/1000  << " " << t8/1000  << " " << t9/1000 << std::endl;
#endif
}


/**
 * Calculates Q*A in place. To construct Q pass identity matrix as input A.
 * @param packed Packed result of tridiagonalization that contains householder vectors that define Q in columns bellow the diagonal. Usually result of a call to `block_householder_tridiag3`.
 * @param[in,out] A On input a matrix to multiply with Q. On output the product Q*A.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
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

/**
 * Calculates Q*A in place. To construct Q pass an appropriate identity matrix as input A.
 * @param packed Packed result of tridiagonalization that contains householder vectors that define Q in columns bellow the diagonal. Usually result of a call to `block_householder_tridiag3`.
 * @param[in,out] A On input a matrix to multiply with Q. On output the product Q*A.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes larger value is optimal.
 */
void block_apply_packed_Q_gpu2(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100) {
  //if input A==Identity, constructs Q
  matrix_gpu A_gpu(A);
  Eigen::MatrixXd scratchSpace(A.rows(), r);
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (int k = (packed.rows() - 3) / r * r; k >= 0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k - 1, actual_r);
    W.col(0) = packed.col(k).tail(W.rows());
    for (size_t j = 1; j < actual_r; j++) {
#ifdef TIME_IT
      start = std::chrono::steady_clock::now();
#endif
      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
#ifdef TIME_IT
      t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      start = std::chrono::steady_clock::now();
#endif
      W.col(j).noalias() = -W.leftCols(j) * scratchSpace.col(0).head(j);
      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
#ifdef TIME_IT
      t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    Eigen::MatrixXd packed_block_transpose_triang = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>();
    matrix_gpu packed_block_transpose_triang_gpu(packed_block_transpose_triang);
    matrix_gpu A_bottom_gpu(A.rows() - k - 1, A.cols());
    A_bottom_gpu.sub_block(A_gpu, k+1, 0, 0, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
    matrix_gpu W_gpu(W);
#ifdef TIME_IT
    t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
    A_bottom_gpu = A_bottom_gpu - W_gpu * (packed_block_transpose_triang_gpu * A_bottom_gpu);
    A_gpu.sub_block(A_bottom_gpu, 0, 0, k+1, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
//    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
//    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
  }
  copy(A,A_gpu);
#ifdef TIME_IT
  std::cout << "block_apply_packed_Q_gpu2 detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

/**
 * Calculates Q*A in place. To construct Q pass an appropriate identity matrix as input A.
 * @param packed Packed result of tridiagonalization that contains householder vectors that define Q in columns bellow the diagonal. Usually result of a call to `block_householder_tridiag3`.
 * @param[in,out] A On input a matrix to multiply with Q. On output the product Q*A.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes larger value is optimal.
 */
void block_apply_packed_Q_gpu3(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100) {
  matrix_gpu A_gpu(A);
  matrix_gpu packed_gpu(packed);
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0;
  auto start = std::chrono::steady_clock::now();
#endif
//  Eigen::MatrixXd scratchSpace(A.rows(), r);
  for (int k = (packed.rows() - 3) / r * r; k >= 0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
//    Eigen::MatrixXd W(packed.rows() - k - 1, actual_r);
    matrix_gpu W_gpu(packed.rows() - k - 1, actual_r);
//    W.col(0) = packed.col(k).tail(W.rows());
    W_gpu.sub_block(packed_gpu, packed.rows() - W_gpu.rows(), k, 0, 0, W_gpu.rows(), 1);
    for (size_t j = 1; j < actual_r; j++) {
//      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
//      W.col(j).noalias() = -W.leftCols(j) * scratchSpace.col(0).head(j);
//      Eigen::VectorXd tmp0 = scratchSpace.col(0).head(j);
//      Eigen::VectorXd tmp1 = W.col(j);
//      Eigen::VectorXd tmp2 = packed.col(j + k).tail(packed.rows() - k - j - 1);
//      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);

      matrix_gpu temp(j,1);
      try {
#ifdef TIME_IT
        start = std::chrono::steady_clock::now();
#endif
        opencl_kernels::eigendecomp_apply_Q1(
                cl::NDRange(j*64), cl::NDRange(64),
                packed_gpu.buffer(), temp.buffer(), packed_gpu.rows(),
                packed_gpu.cols(), k + j + 1, k, k + j);
#ifdef TIME_IT
        t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
        opencl_kernels::eigendecomp_apply_Q2(
                cl::NDRange((W_gpu.rows()+63)/64*64), cl::NDRange(64),
                W_gpu.buffer(), temp.buffer(), packed_gpu.buffer(),
                W_gpu.rows(), W_gpu.cols(), j, k);
#ifdef TIME_IT
        t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
      } catch (cl::Error& e) {
        check_opencl_error("block_apply_packed_Q_gpu3", e);
      }

//      matrix_gpu packed_block(packed.rows() - k - j - 1, j);
//      packed_block.sub_block(packed_gpu, k + j + 1, k, 0, 0, packed.rows() - k - j - 1, j);
//      matrix_gpu packed_col(packed.rows() - k - j - 1, 1);
//      packed_col.sub_block(packed_gpu, k + j + 1, j + k, 0, 0, packed.rows() - k - j - 1, 1);
//      matrix_gpu W_left(W_gpu.rows(), j);
//      W_left.sub_block(W_gpu,0,0,0,0,W_gpu.rows(),j);
//
//      //matrix_gpu W_col = -1 * (W_left * (transpose(packed_block) * packed_col));
//      matrix_gpu W_col = -1 * (W_left * temp);
//      matrix_gpu W_tmp(W_gpu.rows() - j, 1);
//      W_tmp.sub_block(W_col, j, 0, 0, 0, W_gpu.rows() - j, 1);
//      W_tmp = W_tmp + packed_col;
//      W_gpu.sub_block(W_col, 0, 0, 0, j, j, 1);
//      W_gpu.sub_block(W_tmp, 0, 0, j, j, W_gpu.rows() - j, 1);

    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    Eigen::MatrixXd packed_block_transpose_triang = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>();
    matrix_gpu packed_block_transpose_triang_gpu(packed_block_transpose_triang);
#ifdef TIME_IT
    t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
    matrix_gpu A_bottom_gpu(A.rows() - k - 1, A.cols());
    A_bottom_gpu.sub_block(A_gpu, k + 1, 0, 0, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
    A_bottom_gpu = A_bottom_gpu - W_gpu * (packed_block_transpose_triang_gpu * A_bottom_gpu);
    A_gpu.sub_block(A_bottom_gpu, 0, 0, k + 1, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
//    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
//    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
  copy(A, A_gpu);
#ifdef TIME_IT
  std::cout << "block_apply_packed_Q_gpu3 detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

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
 * all elements except diagonal and subdiagonal are 0),T - shift * I = L * D * L^T. Also calculates element growth of D: sum(abs(D)) / abs(sum(D)).
 * @param diag Diagonal of T
 * @param subdiag Subdiagonal of T.
 * @param shift Shift.
 * @param[out] l Subdiagonal of L.
 * @param[out] d_plus Diagonal of D.
 * @return Element growth.
 */
double get_ldl(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, double shift, Eigen::VectorXd& l, Eigen::VectorXd& d_plus) {
  d_plus[0] = diag[0] - shift;
  double element_growth = abs(d_plus[0]);
  double element_growth_denominator = d_plus[0];
  for (int i = 0; i < subdiag.size(); i++) {
    l[i] = subdiag[i] / d_plus[i];
    d_plus[i] *= get_random_perturbation_multiplier();
    d_plus[i + 1] = diag[i + 1] - shift - l[i] * subdiag[i];
    l[i] *= get_random_perturbation_multiplier();
    element_growth += abs(d_plus[i+1]);
    element_growth_denominator += d_plus[i+1];
  }
  d_plus[subdiag.size()] *= get_random_perturbation_multiplier();
  return element_growth / abs(element_growth_denominator);
}

/**
 * Shifts a LDL decomposition. The algorithm is sometimes called stationary quotients-differences with shifts (stqds).
 * D and D+ are diagonal, L and L+ are lower unit triangular (diagonal elements are 1,
 * all elements except diagonal and subdiagonal are 0). L * D * L^T - shift * I = L+ * D * L+^T.
 * Also calculates element growth of D+: sum(abs(D+)) / abs(sum(D+)).
 * @param d Diagonal of D.
 * @param l Subdiagonal of L.
 * @param shift Shift.
 * @param[out] l_plus Subdiagonal of L+.
 * @param[out] d_plus Diagonal of D+.
 * @return Element growth.
 */
double get_perturbed_shifted_ldl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& d_plus) {
  int n = l.size();
  double s = -shift;
  double element_growth = 0;
  double element_growth_denominator = 0;
  for (int i = 0; i < n; i++) {
    double di_perturbed = d[i];
    double li_perturbed = l[i];
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

/**
 * Calculates shifted LDL and UDU factorizations. Combined with twist index they form twisted factorization for calculation
 * of an eigenvector corresponding to eigenvalue that is equal to the shift. Tha algorithm is sometimes called diferential twisted quotient-differences with shifts (dtwqds).
 * L * D * L^T - shift * I = L+ * D+ * L+^T = U- * D- * U-^T
 * D, D+ and D- are diagonal, L and L+ are lower unit triangular (diagonal elements are 1, all elements except diagonal and subdiagonal are 0),
 * U- is upper unit triangular (diagonal elements are 1, all elements except diagonal and superdiagonal are 0)
 * @param d Diagonal of D.
 * @param l Subdiagonal of L.
 * @param shift Shift.
 * @param[out] l_plus Subdiagonal of L+.
 * @param[out] u_minus Superdiagonal of U-.
 * @return Twist index.
 */
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

/**
 * Calculates Sturm count of a LDL decomposition of a tridiagonal matrix - number of eigenvalues larger or equal to shift.
 * Uses stqds - calculation of shifted LDL decomposition algorithm and counts number of positive elements in D.
 * @param d Diagonal of D.
 * @param l Subdiagonal of L.
 * @param shift Shift.
 * @return Sturm count.
 */
int getSturmCountLdl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift) {
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

/**
 * Refines bounds on the i-th largest eigenvalue of LDL decomposition of a matrix using bisection.
 * @param d Diagonal of D.
 * @param l Subdiagonal of L.
 * @param low[in,out] Low bound on the eigenvalue.
 * @param high[in,out] High bound on the eigenvalue.
 * @param i i-th eigenvalue
 */
void eigenvalBisectRefine(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double& low, double& high, int i) {
  int n = d.size();
  double eps = 3e-16;
  while (abs((high - low) / (high + low)) > eps && abs(high - low) > std::numeric_limits<double>::min()) { // second term is for the case where the eigenvalue is 0 and division yields NaN
    double mid = (high + low) * 0.5;
    if (getSturmCountLdl(d, l, mid) > i) {
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
  min_eigval = diag[0] - abs(subdiag[0]);
  max_eigval = diag[0] + abs(subdiag[0]);
  for (int i = 1; i < n - 1; i++) {
    min_eigval = std::min(min_eigval, diag[i] - abs(subdiag[i]) - abs(subdiag[i - 1]));
    max_eigval = std::max(max_eigval, diag[i] + abs(subdiag[i]) + abs(subdiag[i - 1]));
  }
  min_eigval = std::min(min_eigval, diag[n - 1] - abs(subdiag[n - 2]));
  max_eigval = std::max(max_eigval, diag[n - 1] + abs(subdiag[n - 2]));
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

/*
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
*/
struct bisectionTask{
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
  while(!tQueue.empty()){
    int n_valid = std::min(BISECT_K, static_cast<int>(tQueue.size()));
    Eigen::Array<double, BISECT_K, 1> shifts;
    bisectionTask t[BISECT_K];
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
      if(counts[i]>=t[i].start+1){
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
      if(counts[i] <= my_end - 1){
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

/*
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
*/

/**
 * Calculates an eigenvector from twisted factorization T - shift * I = L+ * D+ * L+^T = U- * D- * U-^T.
 * @param l_plus Subdiagonal of the L+.
 * @param u_minus Superdiagonal of the U-.
 * @param subdiag Subdiagonal of T
 * @param i At which column of `eigenvecs` to store resulting vector.
 * @param twist_idx Twist index.
 * @param[out] eigenvecs Matrix in which to sstore resulting vector.
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
  vec /= vec.norm();
}

/**
 * Finds good shift and shifts a LDL decomposition so as to keep element growth low. L * D * L^T - shift * I = L2 * D2 * L2^T.
 * @param l Subdiagonal of L.
 * @param d Diagonal of D.
 * @param low Low bound on wanted shift.
 * @param high High bound on wanted shift.
 * @param max_ele_growth Maximum desired element growth. If no better options are found it might be exceeded.
 * @param max_shift Maximal difference of shhift from wanted bounds.
 * @param l2 Subdiagonal of L2.
 * @param d2 Diagonal of D2.
 * @param shift Shift.
 * @param min_element_growth Element growth achieved with resulting shift.
 */
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
    //cout << " element growth: " << element_growth << " at " << sh << endl;
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
//  cout << "\t\t" << " element growth: " << min_element_growth << " at " << shift << endl;
}
/*
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
  Eigen::VectorXd l2(n - 1), d2(n), l3(n - 1), d3(n), l_plus(n - 1), u_minus(n - 1);

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
//    if (i > cluster_end) {
//      constructAndRefineCluster(l, d, i, high, low, cluster_start, cluster_end);
//      while (true) {
//        // perturb the representation if any eigenvalues in the cluster are too close
//        if (isClusterSeparable(high, low, cluster_start, cluster_end)) {
//          break;
//        }
//        cout << "perturbing at i = " << i << " cluster: " << cluster_start << ", " << cluster_end << endl;
//        perturbRepresentation(l0, d0, l, d);
//        for (int j = cluster_start; j <= cluster_end; j++) {
//          low[j] = low[j] * (1 - copysign(perturbation_range * n, low[j]));
//          high[j] = high[j] * (1 + copysign(perturbation_range * n, high[j]));
//          eigenvalBisectRefine(d, l, low[j], high[j], j);
//        }
//      }
//      if (cluster_end != n - 1) {
//        int next = cluster_end + 1;
//        low[next] = low[next] * (1 - copysign(perturbation_range * n, low[next]));
//        high[next] = high[next] * (1 + copysign(perturbation_range * n, high[next]));
//        eigenvalBisectRefine(d, l, low[next], high[next], next);
//      }
//      min_element_growth = std::numeric_limits<double>::infinity();
//    }
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
 */
/*
Eigen::Array<int, BISECT_K, 1> getSturmCountLdlVec(const Eigen::VectorXd& d, const Eigen::VectorXd& l, Eigen::Array<double, BISECT_K, 1> shifts, int n_valid){
  int n = l.size();
  Eigen::Array<double, BISECT_K, 1> s;
  s.head(n_valid) = -shifts.head(n_valid);
  Eigen::Array<double, BISECT_K, 1> d_plus, l_plus;
  Eigen::Array<bool, BISECT_K, 1> cond;
  Eigen::Array<int, BISECT_K, 1> counts;
  counts.head(n_valid) = Eigen::Array<int, BISECT_K, 1>(n_valid,0);
  for (int i = 0; i < n; i++) {
    d_plus.head(n_valid) = s.head(n_valid) + d[i];
    counts.head(n_valid) += (d_plus.head(n_valid) >= 0).cast<int>();
    cond.head(n_valid) = d_plus.head(n_valid).unaryExpr([](double x){return is_inf(x);}).cast<bool>() &&
                              s.head(n_valid).unaryExpr([](double x){return is_inf(x);}).cast<bool>();
    s.head(n_valid) = l[i] * l[i] * s.head(n_valid) * (d[i] / d_plus.head(n_valid)) - shifts.head(n_valid);
    for(int j=0;j<n_valid;j++){
      if(cond[j]){
        s[j] = l[i] * l[i] * d[i] - shifts[j];
      }
    }
  }
  return counts;
}*/

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
double findInitialShift(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, Eigen::VectorXd& l0, Eigen::VectorXd& d0, double min_eigval, double max_eigval, double max_ele_growth){
  double shift = 0;
  if(min_eigval>0){
    shift=0.9*min_eigval;
  }
  else if(max_eigval<=0){
    shift=0.9*max_eigval;
  }
  double element_growth = get_ldl(diag, subdiag, shift, l0, d0);
  if(element_growth<max_ele_growth){
    return shift;
  }
  double plus=(max_eigval-min_eigval)*1e-15;
  while(!(element_growth<max_ele_growth)){ //if condition is fliped it would be wrong for the case where element_growth is nan
    plus *= -2;
    element_growth = get_ldl(diag, subdiag, shift + plus, l0, d0);
  }
  return shift + plus;
}

struct mrrrTask{
    int start, end;
    double shift; //total shift, not just last one
    Eigen::VectorXd l, d;
    int level;
};

 /**
  * Calculates eigenvalues and eigenvectors of a (preferrably irreducible) tridiagonal matrix T using MRRR algorithm.
  * @param diag Diagonal of of T.
  * @param subdiag Subdiagonal of T.
  * @param eigenvals[out] Eigenvlues.
  * @param eigenvecs[out] Eigenvectors.
  * @param min_rel_sep Minimal relative separation of eigenvalues before computing eigenvectors.
  * @param max_ele_growth Maximal desired element growth of LDL decompositions.
  */
void mrrr2(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, Eigen::Ref<Eigen::VectorXd> eigenvals,
           Eigen::Ref<Eigen::MatrixXd> eigenvecs, double min_rel_sep = 1e-1, double max_ele_growth = 2) {
  double shift_error = 1e-14;
  int n = diag.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  getGresgorin(diag, subdiag, min_eigval, max_eigval);
  Eigen::VectorXd l(n - 1), d(n), l0(n - 1), d0(n);
  double shift0 = findInitialShift(diag, subdiag,l0, d0,min_eigval, max_eigval,max_ele_growth);
//  cout << "init ele growth " << d0.array().abs().sum() / abs(d0.array().sum()) << endl;
//  cout << "shift0 " << shift0 << endl;
  for (int i = 0; i < n; i++) {
    if (i != n - 1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i] = d0[i] * get_random_perturbation_multiplier();
  }
  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();

  auto start = std::chrono::steady_clock::now();
  eigenvalsBisect4(diag, subdiagSquared, min_eigval, max_eigval, low, high);
  std::cout << "eigenvals: "
      << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
      << "ms" << endl;
  eigenvals = (high + low) * 0.5;
  low.array() -= shift0;
  high.array() -= shift0;
  start = std::chrono::steady_clock::now();
  for (int i = 0; i < n; i++) {
    low[i] = low[i] * (1 - copysign(perturbation_range * n, low[i]));
    high[i] = high[i] * (1 + copysign(perturbation_range * n, high[i]));
    eigenvalBisectRefine(d, l, low[i], high[i], i);
  }
   std::cout << "refine0 (shift0 = " << shift0 << "(" << min_eigval << " - " << max_eigval << ")): "
             << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
             << "ms" << endl;
  int t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;
  std::queue<mrrrTask> blockQueue;
  blockQueue.push(mrrrTask{0, n, shift0, std::move(l), std::move(d), 0});
  l.resize(n-1); //after move out
  d.resize(n);
  while (!blockQueue.empty()) {
    mrrrTask block = blockQueue.front();
    blockQueue.pop();
//    cout << "at level "<< block.level << " (" << block.start << ", " << block.end << ")" << endl;
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
      if(cluster_end-i>0){//cluster
        double max_shift = (high[i] - low[cluster_end])*10;
        double currentShift, min_ele_growth;
        findShift(block.l, block.d, low[cluster_end], high[i],max_ele_growth, max_shift, l, d, currentShift, min_ele_growth);
        for(int j=i;j<=cluster_end;j++){
//          if(j >= getSturmCountLdl(block.d, block.l, low[j]) || j < getSturmCountLdl(block.d, block.l, high[j])){
//            cout << "PRE-REFINE ERROR!!!!!" << endl;
//          }
          low[j] = low[j] * (1 - copysign(shift_error, low[j])) - currentShift;
          high[j] = high[j] * (1 + copysign(shift_error, high[j])) - currentShift;
//          if(j >= getSturmCountLdl(d, l, low[j]) || j < getSturmCountLdl(d, l, high[j])){
//            cout << "PRE-REFINE ERROR!!!!!" << endl;
//          }
          eigenvalBisectRefine(d, l, low[j], high[j], j);
//          if(j >= getSturmCountLdl(d, l, low[j]) || j < getSturmCountLdl(d, l, high[j])){
//            cout << "REFINE ERROR!!!!!" << endl;
//          }
        }
//        if(cluster_end+1 == getSturmCountLdl(d, l, low[cluster_end]) && 1 == getSturmCountLdl(d, l, high[i])){
//          cout << "BLOCK ERROR!!!!!" << endl;
//        }
        blockQueue.push(mrrrTask{i,cluster_end+1,block.shift + currentShift, std::move(l), std::move(d), block.level+1});
        l.resize(n-1); //after move out
        d.resize(n);

        i=cluster_end;
      }
      else{ //isolated eigenvalue
//        cout << "\t\t\t\t\t\t\t\tgetting eigenvector " << i << " (eigenvalue " << eigenvals[i] << ")" << endl;
        double low_t=low[i], high_t=high[i];
//        if(i+1 != getSturmCountLdl(block.d, block.l, low[i]) || i != getSturmCountLdl(block.d, block.l, high[i])){
//          cout << "SINGLE ERROR!!!!!" << endl;
//        }
        int twist_idx;
        double low_gap = i == block.start ? std::numeric_limits<double>::infinity() : low[i - 1] - high[i];
        double high_gap = i == block.end - 1 ? std::numeric_limits<double>::infinity() : low[i] - high[i + 1];
        double min_gap = std::min(low_gap, high_gap);
        if (abs(min_gap / ((high[i] + low[i]) * 0.5)) > min_rel_sep) {
          start = std::chrono::steady_clock::now();
          twist_idx = get_twisted_factorization(block.d, block.l, (low[i]+high[i])*0.5, l_plus, u_minus);
          t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          cout << "\t\t" << i << " UNSHIFTED gap: " << min_gap / ((high[i] + low[i]) * 0.5) << endl;
        }
        else if (abs(min_gap / ((high[i] + low[i]) * 0.5 - shift)) > min_rel_sep && min_element_growth < max_ele_growth) {
          low[i] = low[i] * (1 - copysign(shift_error, low[i])) - shift;
          high[i] = high[i] * (1 + copysign(shift_error, high[i])) - shift;
//          if(i >= getSturmCountLdl(d2, l2, low[i]) || i < getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          start = std::chrono::steady_clock::now();
          eigenvalBisectRefine(d2, l2, low[i], high[i], i);
          t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          if(i+1 != getSturmCountLdl(d2, l2, low[i]) || i != getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          double shifted_eigenval = (low[i] + high[i]) * 0.5;
          start = std::chrono::steady_clock::now();
          twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
          t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          cout << "\t\t" << i << " prev shift gap: " << min_gap / ((high[i] + low[i]) * 0.5 - shift) << endl;
        }
        else {
          double max_shift = min_gap / min_rel_sep;
//          if(max_shift<=0){
//            cout << "ERROR!!!!!" << endl;
//          }
          start = std::chrono::steady_clock::now();
          findShift(block.l, block.d, low[i], high[i], max_ele_growth, max_shift, l2, d2, shift, min_element_growth);
          t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();

          //shift and refine eigenvalue
          low[i] = low[i] * (1 - copysign(shift_error, low[i])) - shift;
          high[i] = high[i] * (1 + copysign(shift_error, high[i])) - shift;
//          if(i >= getSturmCountLdl(d2, l2, low[i]) || i < getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          start = std::chrono::steady_clock::now();
          eigenvalBisectRefine(d2, l2, low[i], high[i], i);
          t5+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          if(i+1 != getSturmCountLdl(d2, l2, low[i]) || i != getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          double shifted_eigenval = (low[i] + high[i]) * 0.5;
//          cout << "\t\t" << i << " shifted gap: " << min_gap / shifted_eigenval << endl;
          start = std::chrono::steady_clock::now();
          twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
          t6+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        }
        //calculate eigenvector
        start = std::chrono::steady_clock::now();
        calculateEigenvector(l_plus, u_minus, subdiag, i, twist_idx, eigenvecs);
        t7+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      }
    }
  }
#ifdef TIME_IT
   std::cout << "mrrr detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << " " << t5/1000 << " " << t6/1000 << " " << t7/1000 << std::endl;
#endif
}

/**
 * Calculates eigenvalues and eigenvectors of a (preferrably irreducible) tridiagonal matrix T using MRRR algorithm.
 * @param diag Diagonal of of T.
 * @param subdiag Subdiagonal of T.
 * @param eigenvals[out] Eigenvlues.
 * @param eigenvecs[out] Eigenvectors.
 * @param min_rel_sep Minimal relative separation of eigenvalues before computing eigenvectors.
 * @param max_ele_growth Maximal desired element growth of LDL decompositions.
 */
void mrrr_gpu(const Eigen::Ref<const Eigen::VectorXd> diag, const Eigen::Ref<const Eigen::VectorXd> subdiag, Eigen::Ref<Eigen::VectorXd> eigenvals,
           Eigen::Ref<Eigen::MatrixXd> eigenvecs, double min_rel_sep = 1e-1, double max_ele_growth = 2) {
  double shift_error = 1e-14;
  int n = diag.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  getGresgorin(diag, subdiag, min_eigval, max_eigval);
  Eigen::VectorXd l(n - 1), d(n), l0(n - 1), d0(n);
  double shift0 = findInitialShift(diag, subdiag,l0, d0,min_eigval, max_eigval,max_ele_growth);
//  cout << "init ele growth " << d0.array().abs().sum() / abs(d0.array().sum()) << endl;
//  cout << "shift0 " << shift0 << endl;
  cout << diag.sum() << endl;
  cout << (diag.array()+shift0).sum() << endl;
  cout << (diag.array()-shift0).sum() << endl;
  for (int i = 0; i < n; i++) {
    if (i != n - 1) {
      l[i] = l0[i] * get_random_perturbation_multiplier();
    }
    d[i] = d0[i] * get_random_perturbation_multiplier();
  }
  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();

  auto start = std::chrono::steady_clock::now();

//  eigenvalsBisect4(diag, subdiagSquared, min_eigval, max_eigval, low, high);
//  std::cout << "eigenvals: "
//            << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//            << "ms" << endl;
//  eigenvals = (high + low) * 0.5;
//  low.array() -= shift0;
//  high.array() -= shift0;
//  start = std::chrono::steady_clock::now();
//  for (int i = 0; i < n; i++) {
//    low[i] = low[i] * (1 - copysign(perturbation_range * n, low[i]));
//    high[i] = high[i] * (1 + copysign(perturbation_range * n, high[i]));
////    low[i]=min_eigval -shift0;
////    high[i]=max_eigval -shift0;
//    eigenvalBisectRefine(d, l, low[i], high[i], i);
//  }

  matrix_gpu l_gpu(l);
  matrix_gpu d_gpu(d);
  matrix_gpu low_gpu(n,1);
  matrix_gpu high_gpu(n,1);
  try{
    opencl_kernels::eigenvals_bisect(
            cl::NDRange(n),
            l_gpu.buffer(), d_gpu.buffer(), low_gpu.buffer(), high_gpu.buffer(),
            min_eigval - shift0, max_eigval - shift0, n);
  }
  catch (cl::Error& e) {
    check_opencl_error("block_apply_packed_Q_gpu3", e);
  }
  copy(low,low_gpu);
  copy(high,high_gpu);
  eigenvals = (high + low).array() * 0.5 + shift0;
  std::cout << "gpu eigenvals (shift0 = " << shift0 << "(" << min_eigval << " - " << max_eigval << ")): "
            << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
            << "ms" << endl;
  int t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0;
  std::queue<mrrrTask> blockQueue;
  blockQueue.push(mrrrTask{0, n, shift0, std::move(l), std::move(d), 0});
  l.resize(n-1); //after move out
  d.resize(n);
  while (!blockQueue.empty()) {
    mrrrTask block = blockQueue.front();
    blockQueue.pop();
//    cout << "at level "<< block.level << " (" << block.start << ", " << block.end << ")" << endl;
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
      if(cluster_end-i>0){//cluster
        double max_shift = (high[i] - low[cluster_end])*10;
        double currentShift, min_ele_growth;
        findShift(block.l, block.d, low[cluster_end], high[i],max_ele_growth, max_shift, l, d, currentShift, min_ele_growth);
        for(int j=i;j<=cluster_end;j++){
//          if(j >= getSturmCountLdl(block.d, block.l, low[j]) || j < getSturmCountLdl(block.d, block.l, high[j])){
//            cout << "PRE-REFINE ERROR!!!!!" << endl;
//          }
          low[j] = low[j] * (1 - copysign(shift_error, low[j])) - currentShift;
          high[j] = high[j] * (1 + copysign(shift_error, high[j])) - currentShift;
//          if(j >= getSturmCountLdl(d, l, low[j]) || j < getSturmCountLdl(d, l, high[j])){
//            cout << "PRE-REFINE ERROR!!!!!" << endl;
//          }
          eigenvalBisectRefine(d, l, low[j], high[j], j);
//          if(j >= getSturmCountLdl(d, l, low[j]) || j < getSturmCountLdl(d, l, high[j])){
//            cout << "REFINE ERROR!!!!!" << endl;
//          }
        }
//        if(cluster_end+1 == getSturmCountLdl(d, l, low[cluster_end]) && 1 == getSturmCountLdl(d, l, high[i])){
//          cout << "BLOCK ERROR!!!!!" << endl;
//        }
        blockQueue.push(mrrrTask{i,cluster_end+1,block.shift + currentShift, std::move(l), std::move(d), block.level+1});
        l.resize(n-1); //after move out
        d.resize(n);

        i=cluster_end;
      }
      else{ //isolated eigenvalue
//        cout << "\t\t\t\t\t\t\t\tgetting eigenvector " << i << " (eigenvalue " << eigenvals[i] << ")" << endl;
        double low_t=low[i], high_t=high[i];
//        if(i+1 != getSturmCountLdl(block.d, block.l, low[i]) || i != getSturmCountLdl(block.d, block.l, high[i])){
//          cout << "SINGLE ERROR!!!!!" << endl;
//        }
        int twist_idx;
        double low_gap = i == block.start ? std::numeric_limits<double>::infinity() : low[i - 1] - high[i];
        double high_gap = i == block.end - 1 ? std::numeric_limits<double>::infinity() : low[i] - high[i + 1];
        double min_gap = std::min(low_gap, high_gap);
        if (abs(min_gap / ((high[i] + low[i]) * 0.5)) > min_rel_sep) {
          start = std::chrono::steady_clock::now();
          twist_idx = get_twisted_factorization(block.d, block.l, (low[i]+high[i])*0.5, l_plus, u_minus);
          t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          cout << "\t\t" << i << " UNSHIFTED gap: " << min_gap / ((high[i] + low[i]) * 0.5) << endl;
        }
        else if (abs(min_gap / ((high[i] + low[i]) * 0.5 - shift)) > min_rel_sep && min_element_growth < max_ele_growth) {
          low[i] = low[i] * (1 - copysign(shift_error, low[i])) - shift;
          high[i] = high[i] * (1 + copysign(shift_error, high[i])) - shift;
//          if(i >= getSturmCountLdl(d2, l2, low[i]) || i < getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          start = std::chrono::steady_clock::now();
          eigenvalBisectRefine(d2, l2, low[i], high[i], i);
          t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          if(i+1 != getSturmCountLdl(d2, l2, low[i]) || i != getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          double shifted_eigenval = (low[i] + high[i]) * 0.5;
          start = std::chrono::steady_clock::now();
          twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
          t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          cout << "\t\t" << i << " prev shift gap: " << min_gap / ((high[i] + low[i]) * 0.5 - shift) << endl;
        }
        else {
          double max_shift = min_gap / min_rel_sep;
//          if(max_shift<=0){
//            cout << "ERROR!!!!!" << endl;
//          }
          start = std::chrono::steady_clock::now();
          findShift(block.l, block.d, low[i], high[i], max_ele_growth, max_shift, l2, d2, shift, min_element_growth);
          t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();

          //shift and refine eigenvalue
          low[i] = low[i] * (1 - copysign(shift_error, low[i])) - shift;
          high[i] = high[i] * (1 + copysign(shift_error, high[i])) - shift;
//          if(i >= getSturmCountLdl(d2, l2, low[i]) || i < getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          start = std::chrono::steady_clock::now();
          eigenvalBisectRefine(d2, l2, low[i], high[i], i);
          t5+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
//          if(i+1 != getSturmCountLdl(d2, l2, low[i]) || i != getSturmCountLdl(d2, l2, high[i])){
//            cout << "SINGLE ERROR!!!!!" << endl;
//          }
          double shifted_eigenval = (low[i] + high[i]) * 0.5;
//          cout << "\t\t" << i << " shifted gap: " << min_gap / shifted_eigenval << endl;
          start = std::chrono::steady_clock::now();
          twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
          t6+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        }
        //calculate eigenvector
        start = std::chrono::steady_clock::now();
        calculateEigenvector(l_plus, u_minus, subdiag, i, twist_idx, eigenvecs);
        t7+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
      }
    }
  }
#ifdef TIME_IT
  std::cout << "mrrr_gpu detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << " " << t5/1000 << " " << t6/1000 << " " << t7/1000 << std::endl;
#endif
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
void reducibleTridiagEigenSolver(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs, double splitThreshold = 1e-12) {
  int n = diag.size();
  eigenvecs.resize(n,n);
  eigenvals.resize(n);
  int last = 0;
  for (int i = 0; i < subdiag.size(); i++) {
    if (abs(subdiag[i] / diag[i]) < splitThreshold && abs(subdiag[i] / diag[i + 1]) < splitThreshold) {
//      cout << "split: " << i << endl;
      eigenvecs.block(last, i + 1, i + 1 - last, n - i - 1) = Eigen::MatrixXd::Constant(i + 1 - last, n - i - 1, 0);
      eigenvecs.block(i + 1, last, n - i - 1, i + 1 - last) = Eigen::MatrixXd::Constant(n - i - 1, i + 1 - last, 0);
      if(last==i){
        eigenvecs(last,last)=1;
        eigenvals[last]=diag[last];
      }
      else {
        mrrr_gpu(diag.segment(last, i + 1 - last),
              subdiag.segment(last, i - last),
              eigenvals.segment(last, i + 1 - last),
              eigenvecs.block(last, last, i + 1 - last, i + 1 - last));
      }

      last = i + 1;
    }
  }
  if(last==n-1){
    eigenvecs(last,last)=1;
    eigenvals[last]=diag[last];
  }
  else {
    mrrr_gpu(diag.segment(last, n - last),
          subdiag.segment(last, subdiag.size() - last),
          eigenvals.segment(last, n - last),
          eigenvecs.block(last, last, n - last, n - last));
  }
}

/**
 * Calculates eigenvalues and eigenvectors of a symmetric matrix.
 * @param A The matrix
 * @param eigenvals[out] Eigenvalues.b
 * @param eigenvecs[out] Eigenvectors.
 */
void symmetricEigenSolver(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs) {
  Eigen::MatrixXd packed;
  auto start = std::chrono::steady_clock::now();
  block_householder_tridiag4(A, packed);

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
  reducibleTridiagEigenSolver(diag, subdiag, eigenvals, eigenvecs);

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
