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


namespace stan {
namespace math {

void p(const Eigen::MatrixXd& a) {
  std::cout << a << std::endl;
}

void s(const Eigen::MatrixXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}

void s(const Eigen::VectorXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}

void p(const Eigen::VectorXd& a) {
  std::cout << a << std::endl;
}

void p(const matrix_gpu& a) {
  Eigen::MatrixXd b(a.rows(), a.cols());
  copy(b, a);
  std::cout << b << std::endl;
}

void householder_tridiag(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  for (size_t k = 0; k < T.rows() - 2; k++) {
    Eigen::VectorXd householder = T.block(k, k, T.rows() - k, 1);
    double q=householder.tail(householder.size()-1).norm();
    double alpha = -copysign(q, householder[0]);
    double rsq=alpha*(alpha-householder[0]);
    householder[0]=0;
    householder[1]-=alpha;
    householder*=1/sqrt(rsq);
    householder *= SQRT_2 / householder.norm();
    double n=householder.norm();

    Eigen::MatrixXd H = Eigen::MatrixXd::Identity(T.rows(), T.rows());
    H.bottomRightCorner(T.rows() - k, T.cols() - k) -= householder * householder.transpose();

    //householder.tail(householder.size()-1) /= 2*r_n;
    Eigen::RowVectorXd product = householder.transpose() * T.bottomRightCorner(T.rows() - k, T.cols() - k).transpose();
    Eigen::MatrixXd partial_update = householder * product;
    T.bottomRightCorner(T.rows() - k, T.cols() - k) +=householder * (product*householder) * householder.transpose() - partial_update - partial_update.transpose();
    Eigen::MatrixXd Q_partial_update = (Q.rightCols(Q.cols() - k) * householder) * householder.transpose();
    Q.rightCols(Q.cols() - k) -= Q_partial_update;
  }
}

void block_householder_tridiag(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q, int r = 60) {
  T = A;
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  //Eigen::ArrayXd mags = A.triangularView<Eigen::Lower>().colwise().norm(); //householder norms could be calculated ahead of time (only below diadonal elements)
  for (size_t k = 0; k < T.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(T.rows() - k - 2)});
    Eigen::MatrixXd V(T.rows() - k, actual_r); // TODO manjši V
    V.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);

    for (size_t j = 0; j < actual_r; j++) {
      Eigen::VectorXd householder = T.block(k + j + 1, k + j, T.rows() - k - j - 1, 1);
      double alpha = copysign(householder.norm(), householder[1]);
      double r_n = sqrt(0.5*alpha*(alpha - householder[0]));
      householder.tail(householder.size()-1) /= 2*r_n;
      householder[0]=(householder[0]-alpha)/(2*r_n);
      if (householder.rows() != 1) {
        householder *= SQRT_2 / householder.norm();
      }

      V.col(j).tail(V.rows() - j - 1) = householder;
      Eigen::MatrixXd partial_update = householder * (T.block(k + j + 1, k + j, A.rows() - k - j - 1, actual_r - j).transpose() * householder).transpose();
      T.block(k + j + 1, k + j, A.rows() - k - j - 1, actual_r - j) -= partial_update;
      T.block(k + j, k + j + 1, actual_r - j, A.rows() - k - j - 1) -= partial_update.transpose();
    }

    Eigen::MatrixXd& Y = V;
    Eigen::MatrixXd W = V;
    for (size_t j = 1; j < actual_r; j++) {
      W.col(j) = V.col(j) - W.leftCols(j) * (Y.leftCols(j).transpose() * V.col(j));
    }
    Eigen::MatrixXd T_block = T.block(k, k + actual_r + 1, T.rows() - k, T.cols() - k - actual_r - 1);
    Eigen::MatrixXd product = W.transpose() * T_block;
    Eigen::MatrixXd Y_partial_update = Y * product;
    T.block(k, k + actual_r + 1, T.rows() - k, T.cols() - k - actual_r - 1) -= Y_partial_update;
    T.block(k + actual_r + 1, k, T.cols() - k - actual_r - 1, T.rows() - k) -= Y_partial_update.transpose();

    Eigen::MatrixXd Q_partial_update = (Q.rightCols(Q.cols() - k) * W) * Y.transpose();
    Q.rightCols(Q.cols() - k) -= Q_partial_update;
    //Q.bottomRows(Q.cols() - k) -= Q_partial_update.transpose();
  }
}

void block_householder_tridiag2(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q, int r = 60) {
  T = A;
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  //Eigen::ArrayXd mags = A.triangularView<Eigen::Lower>().colwise().norm(); //householder norms could be calculated ahead of time (only below diadonal elements)
  for (size_t k = 0; k < std::min(T.rows() - 1, T.cols()); k += r) {
    int actual_r = std::min({r, static_cast<int>(T.rows() - k - 1)});
    Eigen::MatrixXd V(T.rows() - k - 1, actual_r); // TODO manjši V
    V.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);

    for (size_t j = 0; j < actual_r; j++) {
      Eigen::VectorXd householder = T.block(k + j + 1, k + j, T.rows() - k - j - 1, 1);
      householder[0] -= copysign(householder.norm(), householder[0]);
      if (householder.rows() != 1) {
        householder *= SQRT_2 / householder.norm();
      }

      V.col(j).tail(V.rows() - j) = householder;
      Eigen::MatrixXd partial_update = householder * (T.block(k + j + 1, k + j, A.rows() - k - j - 1, actual_r - j).transpose() * householder).transpose();
      T.block(k + j + 1, k + j, A.rows() - k - j - 1, actual_r - j) -= partial_update;
      T.block(k + j, k + j + 1, actual_r - j, A.rows() - k - j - 1) -= partial_update.transpose();
    }

    Eigen::MatrixXd& Y = V;
    Eigen::MatrixXd W = V;
    for (size_t j = 1; j < actual_r; j++) {
      W.col(j) = V.col(j) - W.leftCols(j) * (Y.leftCols(j).transpose() * V.col(j));
    }
    Eigen::MatrixXd Y_partial_update = Y * (W.transpose() * T.block(k + 1, k + actual_r + 1, T.rows() - k - 1, T.cols() - k - actual_r - 1));
    T.block(k + 1, k + actual_r + 1, T.rows() - k - 1, T.cols() - k - actual_r - 1) -= Y_partial_update;
    T.block(k + actual_r + 1, k + 1, T.cols() - k - actual_r - 1, T.rows() - k - 1) -= Y_partial_update.transpose();

    Eigen::MatrixXd Q_partial_update = (Q.block(0, 1, Q.rows(), Q.cols() - k - 1) * W) * Y.transpose();
    Q.block(0, 1, Q.rows(), Q.cols() - k - 1) -= Q_partial_update;
    Q.block(1, 0, Q.cols() - k - 1, Q.rows()) -= Q_partial_update.transpose();
  }
}

}
}

#endif //STAN_OPENCL
#endif //STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP
