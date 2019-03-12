#ifndef STAN_MATH_PRIM_MAT_FUN_OPENCL_TRIDIAGONALIZATION_HPP
#define STAN_MATH_PRIM_MAT_FUN_OPENCL_TRIDIAGONALIZATION_HPP

#ifdef STAN_OPENCL

#include <Eigen/QR>
#include <iostream>

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/opencl/subtract.hpp>
#include <stan/math/opencl/add.hpp>
#include <stan/math/opencl/transpose.hpp>

#include <stan/math/opencl/kernels/tridiagonalization.hpp>

using namespace std;

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

void p(const matrix_cl& a) {
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
void block_householder_tridiag_cl(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
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
    matrix_cl U_gpu(U.bottomRows(U.rows() - actual_r + 1).eval());
    matrix_cl V_T_gpu(V.bottomRows(V.rows() - actual_r + 1).transpose().eval());
    matrix_cl partial_update_gpu = U_gpu * V_T_gpu;
    Eigen::MatrixXd partial_update(partial_update_gpu.rows(), partial_update_gpu.cols());
    copy(partial_update, partial_update_gpu);
    packed.block(k + actual_r, k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -= partial_update + partial_update.transpose();
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
  packed(packed.rows() - 2, packed.cols() - 1) = packed(packed.rows() - 1, packed.cols() - 2);
#ifdef TIME_IT
  std::cout << "block_householder_tridiag_cl detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

/**
 * Tridiagonalize a symmetric matrix using block Housholder algorithm. A = Q * T * Q^T, where T is tridiagonal and Q is orthonormal.
 * @param A Input matrix
 * @param[out] packed Packed form of the tridiagonal matrix. Elements of the resulting symmetric tridiagonal matrix T are in the diagonal and first superdiagonal.
 * Columns bellow diagonal contain householder vectors that can be used to construct orthogonal matrix Q.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes a larger value is optimal.
 */
void block_householder_tridiag_cl2(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed, int r = 60) {
  matrix_cl packed_gpu(A);
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0, t5=0, t6=0, t7=0, t8=0, t9=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < A.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(A.rows() - k - 2)});

#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    matrix_cl V_gpu(A.rows() - k - 1, actual_r+1);
    V_gpu.zeros<TriangularViewCL::Upper>();
#ifdef TIME_IT
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif

    for (size_t j = 0; j < actual_r; j++) {
#ifdef TIME_IT
      start = std::chrono::steady_clock::now();
#endif
      matrix_cl Uu(j,1), Vu(j,1), q_gpu(1,1);
      try{
//        opencl_kernels::tridiagonalization_householder_v1(
//                cl::NDRange(128), cl::NDRange(128),
//                packed_gpu.buffer(), V_gpu.buffer(), q_gpu.buffer(), Uu.buffer(), Vu.buffer(),
//                packed_gpu.rows(), V_gpu.rows(), j, k);
        opencl_kernels::tridiagonalization_householder(
                cl::NDRange(1024), cl::NDRange(1024),
                packed_gpu.buffer(), V_gpu.buffer(), q_gpu.buffer(),
                packed_gpu.rows(), V_gpu.rows(), j, k);
#ifdef TIME_IT
        t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
        if(j!=0) {
          opencl_kernels::tridiagonalization_v1(
                  cl::NDRange(64 * j), cl::NDRange(64),
                  packed_gpu.buffer(), V_gpu.buffer(), Uu.buffer(), Vu.buffer(),
                  packed_gpu.rows(), V_gpu.rows(), k);
        }
#ifdef TIME_IT
        t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
//        opencl_kernels::tridiagonalization_v2(
//                cl::NDRange((A.rows() - k - j - 1)*64),cl::NDRange(64),
//                packed_gpu.buffer(), V_gpu.buffer(), Uu.buffer(), Vu.buffer(),
//                packed_gpu.rows(), V_gpu.rows(), k, j);
//        opencl_kernels::tridiagonalization_v2_2(
//                cl::NDRange((A.rows() - k - j - 1 + 31)/32*32),cl::NDRange(32),
//                packed_gpu.buffer(), V_gpu.buffer(), Uu.buffer(), Vu.buffer(),
//                packed_gpu.rows(), V_gpu.rows(), k, j);
        opencl_kernels::tridiagonalization_v2_3(
                cl::NDRange((A.rows() - k - j - 1 + 63)/64*64),cl::NDRange(64),
                packed_gpu.buffer(), V_gpu.buffer(), Uu.buffer(), Vu.buffer(),
                packed_gpu.rows(), V_gpu.rows(), k, j);
#ifdef TIME_IT
        t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
//        opencl_kernels::tridiagonalization_v3(
//                cl::NDRange(128),cl::NDRange(128),
//                packed_gpu.buffer(), V_gpu.buffer(), q_gpu.buffer(),
//                packed_gpu.rows(), V_gpu.rows(), k, j);
        opencl_kernels::tridiagonalization_v3_2(
                cl::NDRange(128),cl::NDRange(128),
                packed_gpu.buffer(), V_gpu.buffer(), q_gpu.buffer(),
                packed_gpu.rows(), V_gpu.rows(), k, j);
      }
      catch (cl::Error& e) {
        check_opencl_error("block_apply_packed_Q_cl3", e);
      }
#ifdef TIME_IT
      t5+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    matrix_cl U_gpu(V_gpu.rows() - actual_r + 1, actual_r);
    U_gpu.sub_block(packed_gpu, k + actual_r, k, 0, 0, A.rows() - k - actual_r, actual_r);
//    matrix_cl V_T_gpu(V.bottomRows(V.rows() - actual_r + 1).transpose().eval());
    matrix_cl Vb_gpu(V_gpu.rows() - actual_r + 1, actual_r);
    Vb_gpu.sub_block(V_gpu, actual_r - 1, 0, 0, 0, V_gpu.rows() - actual_r + 1, actual_r);
#ifdef TIME_IT
    t6+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
    matrix_cl partial_update_gpu = U_gpu * transpose(Vb_gpu);
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
      check_opencl_error("block_apply_packed_Q_cl3", e);
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
  std::cout << "block_householder_tridiag_cl2 detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << " " << t5/1000  << " " << t6/1000  << " " << t7/1000  << " " << t8/1000  << " " << t9/1000 << std::endl;
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
void block_apply_packed_Q_cl2(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 200) {
  //if input A==Identity, constructs Q
  auto start = std::chrono::steady_clock::now();
  matrix_cl A_gpu(A);
//  Eigen::MatrixXd packed_transpose_triang = packed.transpose().triangularView<Eigen::StrictlyUpper>();
//  matrix_cl packed_transpose_triang_gpu(packed_transpose_triang);
  Eigen::MatrixXd scratchSpace(A.rows(), r);
#ifdef TIME_IT
  int t0=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
  int t1=0, t2=0, t3=0, t4=0;
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
    matrix_cl packed_block_transpose_triang_gpu(packed_block_transpose_triang);
//    matrix_cl packed_block_transpose_triang_gpu(actual_r, packed.rows() - k - 1);
//    packed_block_transpose_triang_gpu.sub_block(packed_transpose_triang_gpu, k, k+1, 0, 0, actual_r, packed.rows() - k - 1);
    matrix_cl A_bottom_gpu(A.rows() - k - 1, A.cols());
    A_bottom_gpu.sub_block(A_gpu, k+1, 0, 0, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
    matrix_cl W_gpu(W);
#ifdef TIME_IT
    t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
    A_bottom_gpu = A_bottom_gpu - W_gpu * tri_rect_multiply<TriangularViewCL::Upper>(packed_block_transpose_triang_gpu, A_bottom_gpu);
    A_gpu.sub_block(A_bottom_gpu, 0, 0, k+1, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
//    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
//    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
  }
  copy(A,A_gpu);
#ifdef TIME_IT
  std::cout << "block_apply_packed_Q_cl2 detailed timings: " << t0/1000 << " " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

/**
 * Calculates Q*A in place. To construct Q pass an appropriate identity matrix as input A.
 * @param packed Packed result of tridiagonalization that contains householder vectors that define Q in columns bellow the diagonal. Usually result of a call to `block_householder_tridiag3`.
 * @param[in,out] A On input a matrix to multiply with Q. On output the product Q*A.
 * @param r Block size. Affects only performance of the algorithm. Optimal value depends on the size of A and cache of the processor. For larger matrices or larger cache sizes larger value is optimal.
 */
void block_apply_packed_Q_cl3(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 120) {
  matrix_cl A_gpu(A);
  matrix_cl packed_gpu(packed);
#ifdef TIME_IT
  int t1=0, t2=0, t3=0, t4=0;
  auto start = std::chrono::steady_clock::now();
#endif
//  Eigen::MatrixXd scratchSpace(A.rows(), r);
  for (int k = (packed.rows() - 3) / r * r; k >= 0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
//    Eigen::MatrixXd W(packed.rows() - k - 1, actual_r);
    matrix_cl W_gpu(packed.rows() - k - 1, actual_r);
//    W.col(0) = packed.col(k).tail(W.rows());
    W_gpu.sub_block(packed_gpu, packed.rows() - W_gpu.rows(), k, 0, 0, W_gpu.rows(), 1);
    for (size_t j = 1; j < actual_r; j++) {
//      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
//      W.col(j).noalias() = -W.leftCols(j) * scratchSpace.col(0).head(j);
//      Eigen::VectorXd tmp0 = scratchSpace.col(0).head(j);
//      Eigen::VectorXd tmp1 = W.col(j);
//      Eigen::VectorXd tmp2 = packed.col(j + k).tail(packed.rows() - k - j - 1);
//      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);

      matrix_cl temp(j,1);
      try {
#ifdef TIME_IT
        start = std::chrono::steady_clock::now();
#endif
        opencl_kernels::tridiagonalization_apply_Q1(
                cl::NDRange(j*64), cl::NDRange(64),
                packed_gpu.buffer(), temp.buffer(), packed_gpu.rows(),
                packed_gpu.cols(), k + j + 1, k, k + j);
#ifdef TIME_IT
        t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
        start = std::chrono::steady_clock::now();
#endif
        opencl_kernels::tridiagonalization_apply_Q2(
                cl::NDRange((W_gpu.rows()+63)/64*64), cl::NDRange(64),
                W_gpu.buffer(), temp.buffer(), packed_gpu.buffer(),
                W_gpu.rows(), W_gpu.cols(), j, k);
#ifdef TIME_IT
        t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
      } catch (cl::Error& e) {
        check_opencl_error("block_apply_packed_Q_cl3", e);
      }

//      matrix_cl packed_block(packed.rows() - k - j - 1, j);
//      packed_block.sub_block(packed_gpu, k + j + 1, k, 0, 0, packed.rows() - k - j - 1, j);
//      matrix_cl packed_col(packed.rows() - k - j - 1, 1);
//      packed_col.sub_block(packed_gpu, k + j + 1, j + k, 0, 0, packed.rows() - k - j - 1, 1);
//      matrix_cl W_left(W_gpu.rows(), j);
//      W_left.sub_block(W_gpu,0,0,0,0,W_gpu.rows(),j);
//
//      //matrix_cl W_col = -1 * (W_left * (transpose(packed_block) * packed_col));
//      matrix_cl W_col = -1 * (W_left * temp);
//      matrix_cl W_tmp(W_gpu.rows() - j, 1);
//      W_tmp.sub_block(W_col, j, 0, 0, 0, W_gpu.rows() - j, 1);
//      W_tmp = W_tmp + packed_col;
//      W_gpu.sub_block(W_col, 0, 0, 0, j, j, 1);
//      W_gpu.sub_block(W_tmp, 0, 0, j, j, W_gpu.rows() - j, 1);

    }
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
//    Eigen::MatrixXd packed_block_transpose = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>();
//    matrix_cl packed_block_transpose_gpu(packed_block_transpose);

    matrix_cl packed_block_gpu(packed.rows() - k - 1, actual_r);
    packed_block_gpu.sub_block(packed_gpu, k + 1, k, 0, 0, packed.rows() - k - 1, actual_r);
    matrix_cl packed_block_transpose_gpu = transpose(packed_block_gpu);
#ifdef TIME_IT
    t3+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
    matrix_cl A_bottom_gpu(A.rows() - k - 1, A.cols());
    A_bottom_gpu.sub_block(A_gpu, k + 1, 0, 0, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
    A_bottom_gpu = A_bottom_gpu - W_gpu * tri_rect_multiply<TriangularViewCL::Upper>(packed_block_transpose_gpu, A_bottom_gpu);
    A_gpu.sub_block(A_bottom_gpu, 0, 0, k + 1, 0, A_bottom_gpu.rows(), A_bottom_gpu.cols());
//    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
//    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
#ifdef TIME_IT
    t4+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
  copy(A, A_gpu);
#ifdef TIME_IT
  std::cout << "block_apply_packed_Q_cl3 detailed timings: " << t1/1000 << " " << t2/1000 << " " << t3/1000 << " " << t4/1000 << std::endl;
#endif
}

}
}

#endif
#endif
