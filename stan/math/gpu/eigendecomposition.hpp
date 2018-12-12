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

void s(const Eigen::RowVectorXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}

void p(const Eigen::VectorXd& a) {
  std::cout << a << std::endl;
}

void p(const Eigen::RowVectorXd& a) {
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
    double q = householder.tail(householder.size() - 1).squaredNorm();
    double alpha = -copysign(sqrt(q), householder[0]);
    q -= householder[1] * householder[1];
    householder[0] = 0;
    householder[1] -= alpha;
    q+=householder[1]*householder[1];
    householder *= SQRT_2 / sqrt(q);

    Eigen::VectorXd& u = householder;
    Eigen::VectorXd y = T.bottomRightCorner(T.rows() - k, T.cols() - k) * householder;
    double cnst = y.transpose() * u;
    Eigen::VectorXd v = y - 0.5 * cnst * u;

    Eigen::MatrixXd update = u * v.transpose();
    T.bottomRightCorner(T.rows() - k, T.cols() - k) -= update + update.transpose();

    Q.rightCols(Q.cols() - k) -= (Q.rightCols(Q.cols() - k) * householder) * householder.transpose();
  }
}

void householder_tridiag2(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  for (size_t k = 0; k < T.rows() - 2; k++) {
    Eigen::VectorXd householder = T.block(k + 1, k, T.rows() - k - 1, 1);
    double q = householder.squaredNorm();
    double alpha = -copysign(sqrt(q), T(k,k));
    q -= householder[0] * householder[0];
    householder[0] -= alpha;
    q+=householder[0]*householder[0];
    householder *= SQRT_2 / sqrt(q);

    Eigen::VectorXd& u = householder;
    Eigen::VectorXd y = T.bottomRightCorner(T.rows() - k, T.cols() - k - 1) * u;
    double cnst = y.tail(u.size()).transpose() * u;
    Eigen::VectorXd& v = y;
    v.tail(u.size()) -= 0.5 * cnst * u;

    Eigen::MatrixXd update = u * v.transpose();
    T.bottomRightCorner(T.rows() - k -1, T.cols() - k) -= update;
    T.bottomRightCorner(T.rows() - k, T.cols() - k -1) -= update.transpose();

    Q.rightCols(Q.cols() - k - 1) -= (Q.rightCols(Q.cols() - k - 1) * householder) * householder.transpose();
  }
}

void block_householder_tridiag(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q, int r = 60) {
  T = A;
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  Eigen::MatrixXd Q2=Q;
  Eigen::MatrixXd T2=T;
  //Eigen::ArrayXd mags = A.triangularView<Eigen::Lower>().colwise().norm(); //householder norms could be calculated ahead of time (only below diagonal elements)
  for (size_t k = 0; k < T.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(T.rows() - k - 2)});
    Eigen::MatrixXd V(T.rows(), actual_r); // TODO manjši V
    V.triangularView<Eigen::StrictlyUpper>() = Eigen::MatrixXd::Constant(V.rows(), V.cols(), 0);

    for (size_t j = 0; j < actual_r; j++) {
      Eigen::VectorXd householder = T.block(k + j, k + j, T.rows() - k - j, 1);
      double q = householder.tail(householder.size() - 1).squaredNorm();
      double alpha = -copysign(sqrt(q), householder[0]);
      q -= householder[1] * householder[1];
      householder[0] = 0;
      householder[1] -= alpha;
      q+=householder[1]*householder[1];
      householder *= SQRT_2 / sqrt(q);

      V.col(j).tail(V.rows() - j - k) = householder;
      Eigen::MatrixXd product = T.bottomRightCorner(T.rows() - k - j, T.cols() - k - j) * householder;
      Eigen::RowVectorXd product3 = householder.transpose() * T.block(k + j, k + j, A.rows() - k - j, actual_r - j);
      Eigen::MatrixXd partial_update = householder * product3;
      Eigen::MatrixXd partial_update_t = (product * householder.head(actual_r - j).transpose());
      Eigen::MatrixXd partial_update_q = householder * (householder.transpose() * product) * householder.head(actual_r - j).transpose();

      cout << "partial_update\n" << partial_update << endl;
      //T.block(k + j, k + j, A.rows() - k - j, actual_r - j) -= partial_update;
      //T.block(k + j, k + j, actual_r - j, A.rows() - k - j) -= partial_update.transpose();
      //T.block(k + j, k + j, A.rows() - k - j, A.rows() - k - j) += householder * (product * householder.tail(actual_r-j)) * householder.transpose();

      T.block(k + j, k + j, A.rows() - k - j, actual_r - j) += partial_update_q - partial_update - partial_update_t;

      size_t k2=k + j;
      Eigen::RowVectorXd product2 = householder.transpose() * T2.bottomRightCorner(T.rows() - k2, T.cols() - k2).transpose();
      Eigen::MatrixXd partial_update2 = householder * product2;
      cout << "partial_update2\n" << partial_update2 << endl;
      cout << "partial_update_t\n" << partial_update_t << endl;
      cout << "partial_update_t2\n" << partial_update2.transpose() << endl;
      cout << "partial_update_q\n" << partial_update_q << endl;
      cout << "partial_update_q2\n" << householder * (product2 * householder) * householder.transpose() << endl;
      T2.bottomRightCorner(T2.rows() - k2, T2.cols() - k2) +=
              householder * (product2 * householder) * householder.transpose() -
              partial_update2 - partial_update2.transpose();
      Q2.rightCols(Q2.cols() - k2) -= (Q2.rightCols(Q2.cols() - k2) * householder) * householder.transpose();

      cout << "reconstruct2: " << (A - Q2 * T2 * Q2.transpose()).array().abs().sum() << endl;
      cout << "DT\n" << T2-T << endl;
      int lkh;
    }

    Eigen::MatrixXd& Y = V;
    Eigen::MatrixXd W = V;
    for (size_t j = 1; j < actual_r; j++) {
      W.col(j) = V.col(j) - W.leftCols(j) * (Y.leftCols(j).transpose() * V.col(j));
    }

    /*Eigen::MatrixXd T_block = T.block(k, k + actual_r, T.rows() - k, T.cols() - k - actual_r);
    Eigen::MatrixXd product = W.transpose() * T_block;
    Eigen::MatrixXd T_partial_update = Y * product;
    T.block(k, k + actual_r, T.rows() - k, T.cols() - k - actual_r) -= T_partial_update;
    T.block(k + actual_r, k, T.cols() - k - actual_r, T.rows() - k) -= T_partial_update.transpose();
    Eigen::MatrixXd leftCols = Y.leftCols(Y.cols() - actual_r);
    Eigen::MatrixXd bottomRows = Y.bottomRows(Y.rows() - actual_r);
    T.block(k , k, T.cols() - k, T.rows() - k) -= Y * (product * bottomRows) * W.transpose();*/

    Eigen::MatrixXd T_block = T;//.block(k, k, T.rows() - k, T.cols() - k);
    Eigen::MatrixXd product = W.transpose() * T_block;
    Eigen::MatrixXd T_partial_update = Y * product;
    cout << "partial_update second\n" << T_partial_update << endl;
    T_partial_update.block(k, k, A.rows() - k, actual_r) =  Eigen::MatrixXd::Constant(A.rows() - k, actual_r, 0);
    cout << "partial_update third\n" << T_partial_update << endl;
    T-=T_partial_update;
    T-=T_partial_update.transpose();
    //T.block(k + actual_r, k + actual_r, T.rows() - k - actual_r, T.cols() - k - actual_r) -= T_partial_update.block(actual_r, actual_r,  T_partial_update.rows() - actual_r, T_partial_update.cols() - actual_r);
    //T.block(k + actual_r, k + actual_r, T.cols() - k - actual_r, T.rows() - k - actual_r) -= T_partial_update.block(actual_r, actual_r,  T_partial_update.rows() - actual_r, T_partial_update.cols() - actual_r).transpose();
    Eigen::MatrixXd leftCols = Y.leftCols(Y.cols() - actual_r);
    Eigen::MatrixXd bottomRows = Y.bottomRows(Y.rows() - actual_r);
    Eigen::MatrixXd product1 = Y * (product * Y) * W.transpose();
    T/*.block(k , k, T.cols() - k, T.rows() - k)*/ -= product1/*.bottomRows(product1.rows()-k)*/;
    W=W.bottomRows(W.rows()-k);
    Y=Y.bottomRows(Y.rows()-k);
    cout << "DT_second\n" << T2-T << endl;

    //T.block(k + actual_r, k + actual_r, T.rows() - k, T.cols() - k - actual_r);
    Eigen::MatrixXd Q_partial_update = (Q.rightCols(Q.cols() - k) * W) * Y.transpose();
    Q.rightCols(Q.cols() - k) -= Q_partial_update;
    cout << "DQ\n" << Q2-Q << endl;
    cout << "reconstruct: " << (A - Q * T * Q.transpose()).array().abs().sum() << endl;
    int a;
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
      Eigen::MatrixXd partial_update = householder *
                                       (T.block(k + j + 1, k + j, A.rows() - k - j - 1, actual_r - j).transpose() *
                                        householder).transpose();
      T.block(k + j + 1, k + j, A.rows() - k - j - 1, actual_r - j) -= partial_update;
      T.block(k + j, k + j + 1, actual_r - j, A.rows() - k - j - 1) -= partial_update.transpose();
    }

    Eigen::MatrixXd& Y = V;
    Eigen::MatrixXd W = V;
    for (size_t j = 1; j < actual_r; j++) {
      W.col(j) = V.col(j) - W.leftCols(j) * (Y.leftCols(j).transpose() * V.col(j));
    }
    Eigen::MatrixXd Y_partial_update =
            Y * (W.transpose() * T.block(k + 1, k + actual_r + 1, T.rows() - k - 1, T.cols() - k - actual_r - 1));
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
