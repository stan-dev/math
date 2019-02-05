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
/*
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
    //Eigen::MatrixXd T2=T;
    //T2.bottomRightCorner(T.rows() - k - 1, T.cols() - k) -= update;
    //T2.bottomRightCorner(T.rows() - k, T.cols() - k - 1) -= update.transpose();
    T.bottomRightCorner(T.rows() - k - 1, T.cols() - k - 1) -= update.rightCols(update.cols() - 1) + update.rightCols(update.cols() - 1).transpose();
    T.block(k + 1, k, T.rows() - k - 1, 1) -= update.leftCols(1);
    T.block(k, k + 1, 1, T.rows() - k - 1) -= update.leftCols(1).transpose();

    Q.rightCols(Q.cols() - k - 1) -= (Q.rightCols(Q.cols() - k - 1) * householder) * householder.transpose();
  }
}
void householder_tridiag3(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
#ifndef SKIP_Q
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
#endif
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
    //Eigen::MatrixXd T2=T;
    //T2.bottomRightCorner(T.rows() - k - 1, T.cols() - k) -= update;
    //T2.bottomRightCorner(T.rows() - k, T.cols() - k - 1) -= update.transpose();
    T.bottomRightCorner(T.rows() - k - 1, T.cols() - k - 1) -= update.rightCols(update.cols() - 1) + update.rightCols(update.cols() - 1).transpose();
    T.block(k + 1, k, T.rows() - k - 1, 1) -= update.col(0);
    T.block(k, k + 1, 1, T.rows() - k - 1) -= update.col(0).transpose();

#ifndef SKIP_Q
    Q.rightCols(Q.cols() - k - 1) -= (Q.rightCols(Q.cols() - k - 1) * householder) * householder.transpose();
#endif
  }
}


void householder_tridiag4(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
#ifndef SKIP_Q
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
#endif
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
    Eigen::VectorXd y = T.bottomRightCorner(T.rows() - k, T.cols() - k).template selfadjointView<Eigen::Lower>() * householder;
    double cnst = y.transpose() * u;
    Eigen::VectorXd v = y - 0.5 * cnst * u;

    Eigen::MatrixXd update = u * v.transpose();
    T.bottomRightCorner(T.rows() - k, T.cols() - k) -= update + update.transpose();

#ifndef SKIP_Q
    Q.rightCols(Q.cols() - k).noalias() -= (Q.rightCols(Q.cols() - k) * householder) * householder.transpose();
#endif
  }
}

void householder_tridiag5(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
#ifndef SKIP_Q
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
#endif
#ifdef TIME_IT
  long long t1=0, t2=0, t3=0;
  auto start = std::chrono::steady_clock::now();
  auto stop = std::chrono::steady_clock::now();
#endif
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
    Eigen::VectorXd y = T.bottomRightCorner(T.rows() - k, T.cols() - k).template selfadjointView<Eigen::Lower>() * householder;
    double cnst = y.transpose() * u;
    Eigen::VectorXd v = y - 0.5 * cnst * u;

    T.bottomRightCorner(T.rows() - k, T.cols() - k).template selfadjointView<Eigen::Lower>().rankUpdate(u, v, -1);

#ifdef TIME_IT
    stop = std::chrono::steady_clock::now();
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    start = stop;
#endif

#ifndef SKIP_Q
    Q.rightCols(Q.cols() - k).noalias() -= (Q.rightCols(Q.cols() - k) * householder) * householder.transpose();
#endif

#ifdef TIME_IT
    stop = std::chrono::steady_clock::now();
    t2+=std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    start = stop;
#endif
  }
#ifdef TIME_IT
  start = std::chrono::steady_clock::now();
#endif
  T=T.template selfadjointView<Eigen::Lower>();
#ifdef TIME_IT
  stop = std::chrono::steady_clock::now();
  cout << "out copy: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
       << "ms" << endl;
  cout << "T: " << t1/1000. << "ms" << endl;
  cout << "Q: " << t2/1000. << "ms" << endl;
#endif
}

void householder_tridiag6(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
#ifndef SKIP_Q
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  Eigen::VectorXd scratchSpace(T.rows());
#endif
#ifdef TIME_IT
  long long t1=0, t2=0, t3=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < T.rows() - 2; k++) {
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    Eigen::VectorXd householder = T.block(k, k, T.rows() - k, 1);
    double q = householder.tail(householder.size() - 1).squaredNorm();
    double alpha = -copysign(sqrt(q), householder[0]);
    q -= householder[1] * householder[1];
    householder[0] = 0;
    householder[1] -= alpha;
    q+=householder[1]*householder[1];
    householder *= SQRT_2 / sqrt(q);

    Eigen::VectorXd& u = householder;
    Eigen::VectorXd y = T.bottomRightCorner(T.rows() - k, T.cols() - k).template selfadjointView<Eigen::Lower>() * householder;
    double cnst = y.transpose() * u;
    Eigen::VectorXd v = y - 0.5 * cnst * u;

    for(size_t i=0; i < u.size(); i++){
      T.col(i + k).tail(u.size() - i) -= u[i] * v.tail(u.size()-i) +
                                         v[i] * u.tail(u.size()-i);
    }

#ifdef TIME_IT
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
#ifndef SKIP_Q
//    Eigen::MatrixXd update = u * v.transpose();
//    Eigen::MatrixXd update2=update+update.transpose();
    //T.bottomRightCorner(T.rows() - k, T.cols() - k).template triangularView<Eigen::Lower>() -= update + update.transpose();
    scratchSpace.noalias() = Q.rightCols(Q.cols() - k - 1) * householder.tail(householder.size() - 1);
    Q.rightCols(Q.cols() - k - 1).noalias() -= scratchSpace * householder.tail(householder.size() - 1).transpose();
#endif
#ifdef TIME_IT
    t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
#ifdef TIME_IT
  start = std::chrono::steady_clock::now();
#endif
  T=T.template selfadjointView<Eigen::Lower>();
#ifdef TIME_IT
  cout << "out copy: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  cout << "T: " << t1/1000. << "ms" << endl;
  cout << "Q: " << t2/1000. << "ms" << endl;
#endif
}

void householder_tridiag7(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
#ifndef SKIP_Q
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  Eigen::VectorXd scratchSpace(T.rows());
#endif
#ifdef TIME_IT
  long long t1=0, t2=0, t3=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < T.rows() - 2; k++) {
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    Eigen::VectorXd householder = T.block(k, k, T.rows() - k, 1);
    double q = householder.tail(householder.size() - 1).squaredNorm();
    double alpha = -copysign(sqrt(q), householder[0]);
    q -= householder[1] * householder[1];
    householder[0] = 0;
    householder[1] -= alpha;
    q+=householder[1]*householder[1];
    q=sqrt(q);
    householder *= SQRT_2 / q;

    Eigen::VectorXd& u = householder;
    Eigen::VectorXd v = T.bottomRightCorner(T.rows() - k, T.cols() - k).template selfadjointView<Eigen::Lower>() * householder;
    //cout << v[0] << v[0] - q / SQRT_2 << " " << v[0] - (q / SQRT_2 + alpha * householder[1]) << " " << v[0] - (q / SQRT_2 - alpha * householder[1]) << endl;
    double cnst = v.transpose() * u;
    v -= 0.5 * cnst * u;
    //p(v);

    for(size_t i=1; i < u.size(); i++){
      T.col(i + k).tail(u.size() - i) -= u[i] * v.tail(u.size()-i) +
                                         v[i] * u.tail(u.size()-i);
    }
    //p(u);
    //cout << T(k + 1, k) << " - " << u[0] * v[1] + v[0] * u[1];
    T(k + 1, k) -= u[0] * v[1] + v[0] * u[1];
    //cout << " = " << T(k + 1, k) << endl;
    T.col(k).tail(u.size()-2) = Eigen::VectorXd::Constant(u.size()-2,0);

#ifdef TIME_IT
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
#ifndef SKIP_Q
//    Eigen::MatrixXd update = u * v.transpose();
//    Eigen::MatrixXd update2=update+update.transpose();
    //T.bottomRightCorner(T.rows() - k, T.cols() - k).template triangularView<Eigen::Lower>() -= update + update.transpose();
    scratchSpace.noalias() = Q.rightCols(Q.cols() - k - 1) * householder.tail(householder.size() - 1);
    Q.rightCols(Q.cols() - k - 1).noalias() -= scratchSpace * householder.tail(householder.size() - 1).transpose();
#endif
#ifdef TIME_IT
    t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
#ifdef TIME_IT
  start = std::chrono::steady_clock::now();
#endif
  T=T.template selfadjointView<Eigen::Lower>();
#ifdef TIME_IT
  cout << "out copy: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  cout << "T: " << t1/1000. << "ms" << endl;
  cout << "Q: " << t2/1000. << "ms" << endl;
#endif
}

void householder_tridiag8(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
#ifndef SKIP_Q
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  Eigen::VectorXd scratchSpace(T.rows());
#endif
#ifdef TIME_IT
  long long t1=0, t2=0, t3=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < T.rows() - 2; k++) {
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    //auto householder = T.block(k, k, T.rows() - k, 1);
    //Eigen::VectorXd
    auto householder = T.col(k).tail(T.rows()-k);
    double q = householder.tail(householder.size() - 1).squaredNorm();
    double alpha = -copysign(sqrt(q), householder[0]);
    q -= householder[1] * householder[1];
    double tmp=householder[0];
    householder[0] = 0;
    householder[1] -= alpha;
    q+=householder[1]*householder[1];
    q=sqrt(q);
    householder *= SQRT_2 / q;

    //cout << "############" << k << endl;
    auto& u = householder;
    //p(static_cast<Eigen::VectorXd>(u));
    //p(T);
    Eigen::VectorXd v = T.bottomRightCorner(T.rows() - k, T.cols() - k).template selfadjointView<Eigen::Lower>() * householder;
    //v[0] = q * SQRT_2 - alpha * householder[1];
    v[0] = q / SQRT_2;// - alpha * householder[1];
    //p(v);
    //cout << v[0] - q / SQRT_2 << endl;
    double cnst = v.transpose() * u;
    v -= 0.5 * cnst * u;
    //p(v);

    for(size_t i=1; i < u.size(); i++){
      T.col(i + k).tail(u.size() - i) -= u[i] * v.tail(u.size()-i) +
                                         v[i] * u.tail(u.size()-i);
    }


#ifdef TIME_IT
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
#ifndef SKIP_Q
//    Eigen::MatrixXd update = u * v.transpose();
//    Eigen::MatrixXd update2=update+update.transpose();
    //T.bottomRightCorner(T.rows() - k, T.cols() - k).template triangularView<Eigen::Lower>() -= update + update.transpose();
    scratchSpace.noalias() = Q.rightCols(Q.cols() - k - 1) * householder.tail(householder.size() - 1);
    Q.rightCols(Q.cols() - k - 1).noalias() -= scratchSpace * householder.tail(householder.size() - 1).transpose();
#endif
#ifdef TIME_IT
    t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    //p(static_cast<Eigen::VectorXd>(u));
    //cout << T(k + 1, k) * q / SQRT_2 + alpha << " - " << (u[0] * v[1] + v[0] * u[1]);
    T(k + 1, k) = T(k + 1, k) * q / SQRT_2 + alpha - v[0] * u[1];
    //cout << " = " << T(k + 1, k) << endl;
    householder[0]=tmp;
    T.col(k).tail(u.size()-2) = Eigen::VectorXd::Constant(u.size()-2,0);
    //T(k + 1, k) -= u[0] * v[1] + v[0] * u[1];
  }
#ifdef TIME_IT
  start = std::chrono::steady_clock::now();
#endif
  T=T.template selfadjointView<Eigen::Lower>();
#ifdef TIME_IT
  cout << "out copy: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  cout << "T: " << t1/1000. << "ms" << endl;
  cout << "Q: " << t2/1000. << "ms" << endl;
#endif
}

void householder_tridiag9(const Eigen::MatrixXd& A, Eigen::MatrixXd& T, Eigen::MatrixXd& Q) {
  //matrix_gpu R_gpu(A);
  //matrix_gpu Q_gpu(Eigen::MatrixXd::Identity(A.rows()));
  T = A;
#ifndef SKIP_Q
  Q = Eigen::MatrixXd::Identity(T.rows(), T.rows());
  Eigen::VectorXd scratchSpace(T.rows());
#endif
#ifdef TIME_IT
  long long t1=0, t2=0, t3=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < T.rows() - 2; k++) {
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    //auto householder = T.block(k, k, T.rows() - k, 1);
    //Eigen::VectorXd
    auto householder = T.col(k).tail(T.rows()-k -1);
    double q = householder.squaredNorm();
    double alpha = -copysign(sqrt(q), T(k,k));
    q -= householder[0] * householder[0];
    householder[0] -= alpha;
    q+=householder[0]*householder[0];
    q=sqrt(q);
    householder *= SQRT_2 / q;

    //cout << "############" << k << endl;
    auto& u = householder;
    //p(static_cast<Eigen::VectorXd>(u));
    //p(T);
    Eigen::VectorXd v(householder.size()+1);
    v.tail(householder.size()).noalias() = T.bottomRightCorner(T.rows() - k - 1, T.cols() - k - 1).template selfadjointView<Eigen::Lower>() * householder;
    //v[0] = q * SQRT_2 - alpha * householder[1];
    v[0] = q / SQRT_2;// - alpha * householder[1];
    //cout << v[0] - q / SQRT_2 << endl;
    //p(v);
    double cnst = v.tail(householder.size()).transpose() * u;
    v.tail(householder.size()) -= 0.5 * cnst * u;
    //p(v);

    for(size_t i=1; i < v.size(); i++){
      T.col(i + k).tail(v.size() - i).noalias() -= u[i-1] * v.tail(v.size()-i) +
                                                   v[i] * u.tail(v.size()-i);
    }


#ifdef TIME_IT
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
    start = std::chrono::steady_clock::now();
#endif
#ifndef SKIP_Q
    //    Eigen::MatrixXd update = u * v.transpose();
//    Eigen::MatrixXd update2=update+update.transpose();
    //T.bottomRightCorner(T.rows() - k, T.cols() - k).template triangularView<Eigen::Lower>() -= update + update.transpose();
    scratchSpace.noalias() = Q.rightCols(Q.cols() - k - 1) * householder;
    Q.rightCols(Q.cols() - k - 1).noalias() -= scratchSpace * householder.transpose();
#endif
#ifdef TIME_IT
    t2+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
    //p(static_cast<Eigen::VectorXd>(u));
    //cout << T(k + 1, k) * q / SQRT_2 + alpha << " - " << (u[0] * v[1] + v[0] * u[1]);
    T(k + 1, k) = T(k + 1, k) * q / SQRT_2 + alpha - v[0] * u[0];
    //cout << " = " << T(k + 1, k) << endl;
    T.col(k).tail(v.size()-2) = Eigen::VectorXd::Constant(v.size()-2,0);
    //T(k + 1, k) -= u[0] * v[1] + v[0] * u[1];
  }
#ifdef TIME_IT
  start = std::chrono::steady_clock::now();
#endif
  T=T.template selfadjointView<Eigen::Lower>();
#ifdef TIME_IT
  cout << "out copy: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  cout << "T: " << t1/1000. << "ms" << endl;
  cout << "Q: " << t2/1000. << "ms" << endl;
#endif
}

*/

/*
void householder_tridiag_packed(const Eigen::MatrixXd& A, Eigen::MatrixXd& packed) {
  packed = A;
#ifdef TIME_IT
  long long t1=0;
  auto start = std::chrono::steady_clock::now();
#endif
  for (size_t k = 0; k < packed.rows() - 2; k++) {
#ifdef TIME_IT
    start = std::chrono::steady_clock::now();
#endif
    //auto householder = packed.block(k, k, packed.rows() - k, 1);
    auto householder = packed.col(k).tail(packed.rows()-k -1);
    double q = householder.squaredNorm();
    double alpha = -copysign(sqrt(q), packed(k,k));
    q -= householder[0] * householder[0];
    householder[0] -= alpha;
    q+=householder[0]*householder[0];
    q=sqrt(q);
    householder *= SQRT_2 / q;

    auto& u = householder;
    Eigen::VectorXd v(householder.size()+1);
    v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - 1, packed.cols() - k - 1).template selfadjointView<Eigen::Lower>() * householder;
    //v[0] = q * SQRT_2 - alpha * householder[1];
    v[0] = q / SQRT_2;
    double cnst = v.tail(householder.size()).transpose() * u;
    v.tail(householder.size()) -= 0.5 * cnst * u;

    for(size_t i=1; i < v.size(); i++){
      packed.col(i + k).tail(v.size() - i) -= u[i-1] * v.tail(v.size()-i) +
                                         v[i] * u.tail(v.size()-i);
    }
    packed(k, k + 1) = packed(k + 1, k) * q / SQRT_2 + alpha - v[0] * u[0];
    //packed.col(k).tail(v.size()-2) = Eigen::VectorXd::Constant(v.size()-2,0);
    //packed(k + 1, k) -= u[0] * v[1] + v[0] * u[1];
#ifdef TIME_IT
    t1+=std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start).count();
#endif
  }
  packed(packed.rows()-2, packed.cols()-1) = packed(packed.rows()-1, packed.cols()-2);
#ifdef TIME_IT
  cout << "packed: " << t1/1000. << "ms" << endl;
#endif
}
*/

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
//      v.tail(householder.size())=(packed.bottomRightCorner(packed.rows()- k - j - 1, packed.cols()- k - j -1) - update - update.transpose())
//                                         .template selfadjointView<Eigen::Lower>() * u;
      v.tail(householder.size()) = packed.bottomRightCorner(packed.rows() - k - j - 1, packed.cols() - k - j -1).selfadjointView<Eigen::Lower>() * u
                                  - U.bottomLeftCorner(u.size(), j) * (V.bottomLeftCorner(u.size(), j).transpose() * u)
                                  - V.bottomLeftCorner(u.size(), j) * (U.bottomLeftCorner(u.size(), j).transpose() * u);
      v[0] = q / SQRT_2;// - alpha * householder[1];
      double cnst = v.tail(householder.size()).transpose() * u;
      v.tail(householder.size()) -= 0.5 * cnst * u;

//
//      for(size_t i=1; i < actual_r - j; i++){//TODO za vsak HH posebaj !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//        packed.col(i + k + j).tail(v.size() - i) -= u[i-1] * v.tail(v.size()-i) +
//                                                v[i] * u.tail(v.size()-i);
//      }
//      p(householder.eval());
//      p(v);
//      p(packed);

      packed(k + j, k + j + 1) = packed(k + j + 1, k + j) * q / SQRT_2 + alpha - v[0] * u[0];
      U.col(j).tail(U.rows() - j) = u;
      V.col(j).tail(V.rows() - j) = v.tail(V.rows() - j);
    }

    Eigen::MatrixXd& Y = U;
//    Eigen::MatrixXd W = U;
//    for (size_t j = 1; j < actual_r; j++) {
//      W.col(j) = U.col(j) - W.leftCols(j) * (U.leftCols(j).transpose() * U.col(j));
//    }
    Eigen::MatrixXd Y_partial_update = U * V.transpose();
    //packed.block(k , k, packed.rows() - k, packed.cols() - k) -= Y_partial_update + Y_partial_update.transpose();
    packed.block(k + actual_r , k + actual_r, packed.rows() - k - actual_r, packed.cols() - k - actual_r).triangularView<Eigen::Lower>() -=
            (Y_partial_update + Y_partial_update.transpose()).bottomRightCorner(Y_partial_update.rows()-actual_r + 1, Y_partial_update.cols()-actual_r + 1);
    //packed.block(k , k, packed.rows() - k, packed.cols() - k).triangularView<Eigen::Lower>() -= U * V.transpose() + V * U.transpose();
    //packed.block(k , k, packed.rows() - k, packed.cols() - k) -= U * V.transpose() + V * U.transpose();
//    for(size_t i=0; i < U.rows(); i+=r) {
//      int actual_r = std::min({r, static_cast<int>(packed.rows() - k - i)});
////      Eigen::MatrixXd T_block = T.block(k + i, k + i, T.rows() - k  - i, actual_r);
////      Eigen::MatrixXd V_block = V.block(i,0,actual_r,V.cols());
////      Eigen::MatrixXd U_block = U.block(i,0,U.rows()-i,U.cols());
////      Eigen::MatrixXd P = U * V.block(i,0,actual_r,V.cols()).transpose();
//      packed.block(k + i, k + i, packed.rows() - k  - i, actual_r) -= U.block(i,0,U.rows()-i,U.cols()) * V.block(i,0,actual_r,V.cols()).transpose() +
//                                                            V.block(i,0,U.rows()-i,U.cols()) * U.block(i,0,actual_r,V.cols()).transpose();
//      //packed.col(k + i).tail(U.rows()-i) = U.
//    }

    //Q.rightCols(Q.cols() - k) -= (Q.rightCols(Q.cols() - k) * W) * Y.transpose();
  }
  packed(packed.rows()-2, packed.cols()-1) = packed(packed.rows()-1, packed.cols()-2);
}

/*
void apply_packed_Q(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A){
  //if input A==Identity, constructs Q
  Eigen::VectorXd scratchSpace(A.rows());
  for (size_t k = 0; k < packed.rows() - 2; k++) {
  //for (int k = packed.rows() - 3; k >=0; k--) {
    auto householder = packed.col(k).tail(packed.rows() - k - 1);
//    p(static_cast<Eigen::VectorXd>(householder));
    scratchSpace.noalias() = A.rightCols(A.cols() - k - 1) * householder;
    A.rightCols(A.cols() - k - 1).noalias() -= scratchSpace * householder.transpose();
//    scratchSpace.head(A.cols() - k - 1).noalias() = A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1) * householder;
//    A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1).noalias() -= scratchSpace.head(A.cols() - k - 1) * householder.transpose();
  }
}

void apply_packed_Q2(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A){
  //if input A==Identity, constructs Q
  A.transposeInPlace();
  Eigen::VectorXd scratchSpace(A.rows());
  //for (size_t k = 0; k < packed.rows() - 2; k++) {
  for (int k = packed.rows() - 3; k >=0; k--) {
    auto householder = packed.col(k).tail(packed.rows() - k - 1);
//    p(static_cast<Eigen::VectorXd>(householder));
//    scratchSpace.noalias() = A.rightCols(A.cols() - k - 1) * householder;
//    A.rightCols(A.cols() - k - 1).noalias() -= scratchSpace * householder.transpose();
    scratchSpace.head(A.cols() - k - 1).noalias() = A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1) * householder;
    A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1).noalias() -= scratchSpace.head(A.cols() - k - 1) * householder.transpose();
  }
  A.transposeInPlace();
}
 */
/*
void apply_packed_Q3(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A){
  //if input A==Identity, constructs Q
  Eigen::RowVectorXd scratchSpace(A.rows());
  //for (size_t k = 0; k < packed.rows() - 2; k++) {
  for (int k = packed.rows() - 3; k >=0; k--) {
    auto householder = packed.col(k).tail(packed.rows() - k - 1);
    scratchSpace.head(A.cols() - k - 1).noalias() = householder.transpose() * A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1);
    A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1).noalias() -= householder * scratchSpace.head(A.cols() - k - 1);
  }
}

void apply_packed_Q4(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A){
  //if input A==Identity, constructs Q
  Eigen::RowVectorXd scratchSpace(A.rows());
  //for (size_t k = 0; k < packed.rows() - 2; k++) {
  for (int k = packed.rows() - 3; k >=0; k--) {
    auto householder = packed.col(k).tail(packed.rows() - k - 1);
    scratchSpace.head(A.cols()).noalias() = householder.transpose() * A.bottomRows(A.cols() - k - 1);
    A.bottomRows(A.cols() - k - 1).noalias() -= householder * scratchSpace.head(A.cols());
  }
}
 */

/*


void apply_packed_Q4(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A){
  //if input A==Identity, constructs Q
  Eigen::VectorXd scratchSpace(A.rows());
  //for (size_t k = 0; k < packed.rows() - 2; k++) {
  for (int k = packed.rows() - 3; k >=0; k--) {
    auto householder = packed.col(k).tail(packed.rows() - k - 1);
    scratchSpace.head(A.cols() - k - 1).noalias() = A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1).transpose() * householder;
    A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1).noalias() -= householder * scratchSpace.head(A.cols() - k - 1).transpose();
  }
}
*/

/*
void block_apply_packed_Q(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100){
  //if input A==Identity, constructs Q
  for (size_t k = 0; k < packed.rows() - 2; k += r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});

    Eigen::MatrixXd U = packed.block(k,k,packed.rows()-k,actual_r).triangularView<Eigen::StrictlyLower>();
    auto& Y = U;
    Eigen::MatrixXd W = U;
    for (size_t j = 1; j < actual_r; j++) {
      W.col(j) = U.col(j) - W.leftCols(j) * (U.leftCols(j).transpose() * U.col(j));
    }
    A.rightCols(A.cols() - k) -= (A.rightCols(A.cols() - k) * W) * Y.transpose();
  }
}
 */

/*
//ne dela
void block_apply_packed_Q2(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100){
  //if input A==Identity, constructs Q
  Eigen::MatrixXd scratchSpace(A.rows(), r);
  for (int k = (packed.rows() - 3)/r*r; k >=0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k - 1, actual_r);
    W.col(0).head(actual_r - 1) = Eigen::VectorXd::Constant(actual_r-1, 0);
    W.col(0).tail(W.rows() - actual_r + 1) = packed.col(k + actual_r - 1).tail(W.rows() - actual_r + 1);
//    for (int j = 1; j < actual_r; j++) {
    for (int j = actual_r-2; j >= 0; j--) {
      scratchSpace.col(0).head(actual_r - j - 1).noalias() = packed.block(k + j + 1, k + j + 1, packed.rows() - k - j - 1, actual_r - j - 1).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1); //TODO za 1 krajsi vec
      W.col(j).noalias() = - W.rightCols(actual_r - j - 1) * scratchSpace.col(0).head(actual_r - j - 1); //TODO .col(j).tail( ... )
      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
    }
    scratchSpace.leftCols(actual_r).noalias() = A.rightCols(A.cols() - k - 1) * W;
    A.rightCols(A.cols() - k - 1).noalias() -= scratchSpace.leftCols(actual_r) * packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>();
  }
}
 */

void block_apply_packed_Q3(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100){
  //if input A==Identity, constructs Q
  Eigen::MatrixXd scratchSpace(A.rows(), r);
  //for (size_t k = 0; k < packed.rows() - 2; k += r) {
  for (int k = (packed.rows() - 3)/r*r; k >=0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k -1, actual_r);
    W.col(0) = packed.col(k).tail(W.rows());
    for (size_t j = 1; j < actual_r; j++) {
      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
      W.col(j).noalias() = - W.leftCols(j) * scratchSpace.col(0).head(j);
      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
    }
//    scratchSpace.leftCols(actual_r).noalias() = A.rightCols(A.cols() - k - 1) * W;
//    A.rightCols(A.cols() - k - 1).noalias() -= scratchSpace.leftCols(actual_r) * packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>();
    scratchSpace.transpose().bottomRows(actual_r).noalias() = packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
    A.bottomRows(A.cols() - k - 1).noalias() -= W * scratchSpace.transpose().bottomRows(actual_r);
    //A.bottomRows(A.rows() - k - 1) -= W * packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>() * A.bottomRows(A.rows() - k - 1);
  }
}

/*
//ne dela
void block_apply_packed_Q3(const Eigen::MatrixXd& packed, Eigen::MatrixXd& A, int r = 100){
  //if input A==Identity, constructs Q
  Eigen::MatrixXd scratchSpace(A.rows(), r);
  //for (size_t k = 0; k < packed.rows() - 2; k += r) {
  for (int k = (packed.rows() - 3)/r*r; k >=0; k -= r) {
    int actual_r = std::min({r, static_cast<int>(packed.rows() - k - 2)});
    Eigen::MatrixXd W(packed.rows() - k -1, actual_r);
//    W.col(actual_r - 1).tail(2) = packed.col(k + actual_r - 1).tail(2);
//    W.col(actual_r - 1).head(W.rows() - 2) = Eigen::VectorXd::Constant(W.rows()-2, 0);
    W.col(0) = packed.col(k).tail(W.rows());

    for (size_t j = 1; j < actual_r; j++) {
//    for (int j = actual_r-2; j >= 0; j--) {
//      scratchSpace.col(0).head(actual_r - j - 1).noalias() = packed.block(k + j + 1, k + j + 1, packed.rows() - k - j - 1, actual_r - j - 1).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
//      W.col(j).noalias() = - W.rightCols(actual_r - j - 1) * scratchSpace.col(0).head(actual_r - j - 1);
//      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
//      W.col(j).head(j) = Eigen::VectorXd::Constant(j, 0);
      scratchSpace.col(0).head(j).noalias() = packed.block(k + j + 1, k, packed.rows() - k - j - 1, j).transpose() * packed.col(j + k).tail(packed.rows() - k - j - 1);
      W.col(j).noalias() = - W.leftCols(j) * scratchSpace.col(0).head(j);
      W.col(j).tail(W.rows() - j) += packed.col(j + k).tail(packed.rows() - k - j - 1);
    }
    Eigen::Map<Eigen::MatrixXd> scratchSpace2(scratchSpace.data(),A.cols() - k - 1,actual_r); //ensure data is contiguous
    scratchSpace2.noalias() = A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1).transpose() * W;
    //Eigen::MatrixXd temp = (scratchSpace2 * packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>()).transpose();
    A.bottomRightCorner(A.cols() - k - 1, A.cols() - k - 1).noalias() -= (scratchSpace2 * packed.block(k + 1, k, packed.rows() - k - 1, actual_r).transpose().triangularView<Eigen::Upper>()).transpose();
  }
}
*/
 //#define DEBUG_PRINT

/*
void dqds_work(Eigen::VectorXd& q, Eigen::VectorXd& e, int start, int stop){
  const double eps=1e-26;

  double shift=0;
  while (start!=stop) {
    //norm = 0;
    double d = q[start] - shift;
    double new_shift = INFTY;
    //unroll while loop and dont check for splits in every iteration
    int short_stop=stop;
    for (int i = start; i < stop; i++) {
      q[i] = d + e[i];
      double temp = q[i + 1] / q[i];
      e[i] *= temp;
      d = d * temp - shift;
      if (d < 0) {
        cout << "ERROR" << endl;
        throw;
      }
      //if(e[i] * e[i] / (q[i] * q[i]) + e[i] * e[i] / (q[i+1] * q[i+1]) < eps){
      if(e[i] * e[i] * (q[i] * q[i] + q[i+1] * q[i+1]) < eps * q[i] * q[i] * q[i+1] * q[i+1]){
        if(i==stop-1){
          short_stop--;
          #ifdef DEBUG_PRINT
          cout << "stop=" << short_stop << endl;
        #endif
        }
        else if (i!=start){
          #ifdef DEBUG_PRINT
          cout << "Split! " << start << " " << i+1 << " " << stop << endl;
          #endif
          dqds_work(q,e,start,i+1);
          //dqds_work(q,e,i+1,stop);
          //return;
          start=i+1;
          continue;
        }
        else{
          start++;
          #ifdef DEBUG_PRINT
          cout << "start=" << start << endl;
        #endif
        }
      }
//      if(i==0){
//        //temp = q[i] + e[i] - sqrt(q[i + 1] * e[i]);
//        temp = q[i] - abs(e[i]);
//      }
//      else {
//        temp = q[i] - abs(e[i]) - abs(e[i-1]);
//        //temp = q[i] + e[i] - sqrt(q[i] * e[i - 1]) - sqrt(q[i + 1] * e[i]);
//      }
//      new_shift=std::min(new_shift,temp);
//      //norm += e[i] * e[i] / (q[i] * q[i]) + e[i] * e[i] / (q[i+1] * q[i+1]);

    }
    q[stop] = d;
    stop=short_stop;
    for(int j=0;j<15;j++) {
      if (start == stop) {
        break;
      }
      d = q[start] - shift;
      for (int i = start; i < stop; i++) {
        q[i] = d + e[i];
        double temp = q[i + 1] / q[i];
        e[i] *= temp;
        d = d * temp - shift;
        if (d < 0) {
          cout << "ERROR" << endl;
          throw;
        }
      }
      q[stop] = d;
      if (e[stop - 1] * e[stop - 1] / (q[stop - 1] * q[stop - 1]) + e[stop - 1] * e[stop - 1] / (q[stop] * q[stop]) < eps) {
        stop--;
#ifdef DEBUG_PRINT
        cout << "stop=" << stop << endl;
#endif
      }
    }
    //shift = std::max(0.,new_shift);
    //cout << "next shift: " << shift << endl;
    //q[stop] -= e[stop-1] - shift;
  }
}

void dqds(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, Eigen::VectorXd& eigenvals){
  int n = subdiag.size();
  //za bidiag
//  Eigen::VectorXd q=diag.array()*diag.array();
//  Eigen::VectorXd e=subdiag.array()*subdiag.array();
  double min_eigval = diag[0] - abs(subdiag[0]); //calc lower bound on eigvals using Gresgorin discs
  for(int i=1;i<n;i++){
    min_eigval = std::min(min_eigval, diag[i] - abs(subdiag[i]) - abs(subdiag[i-1]));
  }
  min_eigval = std::min(min_eigval, diag[n] - abs(subdiag[n-1]));
  cout << "min eigval:" << min_eigval << endl;
  Eigen::VectorXd q(n+1);
  Eigen::VectorXd e(n);
  q[0] = diag[0] - min_eigval;
  for(int i=0; i < n; i++){//TODO is this LDL?
    e[i] = (subdiag[i] / q[i]) * subdiag[i];
    q[i+1] = diag[i+1] - e[i] - min_eigval;
  }

  //double norm=2*eps*n;
  dqds_work(q,e,0,n);
  cout << "end" << endl;
  eigenvals = q.array() + min_eigval;
}
 */


void get_ldl(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, double shift, Eigen::VectorXd& l, Eigen::VectorXd& d_plus){
  d_plus[0] = diag[0] - shift;
  for(int i = 0; i < subdiag.size(); i++){
    l[i] = subdiag[i] / d_plus[i];
    d_plus[i+1] = diag[i+1] - shift - l[i] * subdiag[i];
  }
}

/*
void get_udu(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag, Eigen::VectorXd& u, Eigen::VectorXd& d_minus){
  int n = subdiag.size();
  d_minus[n] = diag[n];
  for(int i = n - 1; i >= 0; i--){
    u[i] = subdiag[i] / d_minus[i+1];
    d_minus[i] = diag[i] - u[i] * subdiag[i];
  }
}
 */

/*
//stationary qds
void get_shifted_ldl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& d_plus, Eigen::VectorXd& s){
  int n = l.size();
  s[0] = -shift;
  for(int i=0; i < n; i++){
    d_plus[i] = s[i] + d[i];
    l_plus[i] = l[i] * d[i] / d_plus[i];
    s[i+1] = l_plus[i] * l[i] * s[i] - shift;
  }
  d_plus[n] = s[n] + d[n];
}
 */

inline double get_random_perturbation_multiplier(){
  static const double perturbation_range = 1e-14;
  static const double rand_norm = perturbation_range / RAND_MAX;
  static const double almost_one = 1 - perturbation_range * 0.5;
  return almost_one + std::rand() * rand_norm;
}

//stationary qds
double get_perturbed_shifted_ldl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& d_plus){
  int n = l.size();
  double s = -shift;
  double element_growth = 0;
  double element_growth_denominator = 0;
  for(int i=0; i < n; i++){
    double di_perturbed = d[i] * get_random_perturbation_multiplier();
    double li_perturbed = l[i] * get_random_perturbation_multiplier();
    d_plus[i] = s + di_perturbed;
    element_growth += abs(d_plus[i]);
    element_growth_denominator += d_plus[i];
    l_plus[i] = li_perturbed * di_perturbed / d_plus[i];
    if(is_inf(d_plus[i]) && is_inf(s)) { // this happens if d_plus[i]==0 -> in next iteration d_plus==inf and s==inf
      s = li_perturbed * li_perturbed * di_perturbed - shift;
    }
    else {
      s = l_plus[i] * li_perturbed * s - shift;
    }
  }
  d_plus[n] = s + d[n] * get_random_perturbation_multiplier();
  element_growth += abs(d_plus[n]);
  return element_growth / abs(element_growth_denominator);
}

/*
//progressive qds
void get_shifted_udu(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& u_minus, Eigen::VectorXd& d_minus, Eigen::VectorXd& p){
  int n = l.size();
  p[n]=d[n]-shift;
  for(int i = n - 1; i >= 0; i--){
    d_minus[i+1] = d[i] * l[i] * l[i] + p[i+1];
    u_minus[i] = l[i] * d[i] / d_minus[i+1];
    p[i] = p[i+1] / d_minus[i+1] * d[i] - shift;
  }
  d_minus[0] = p[0];
}
 */

//dtwqds
int get_twisted_factorization(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift, Eigen::VectorXd& l_plus, Eigen::VectorXd& u_minus){
  int n = l.size();
  //calculate shifted ldl
  Eigen::VectorXd s(n+1);
  s[0] = -shift;
  for(int i=0; i < n; i++){
    double d_plus = s[i] + d[i];
    l_plus[i] = l[i] * d[i] / d_plus;
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
    /*if(is_nan(s[i+1])){
      if(is_inf(d_plus)){ // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
        s[i+1] = l[i] * l[i] * d[i] - shift;
      }
      else if(d_plus==0){
        if(l[i]==0){
          l_plus[i]=d[i];
        }
        else if(d[i]==0){
          l_plus[i]=l[i];
        }
        s[i+1] = l_plus[i] * l[i] * s[i] - shift;
      }
    }*/
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

    /*if(d_minus != 0) {
      p = p * t - shift; //inf                                                              NaN
    }
    else{ // d_minus==0, (t==inf/nan, u_minus[i]==inf/nan, gamma==inf/nan)
      if(abs(d[i]) < zero_threshold){ //t==nan, u_minus[i]==nan, gamma==nan
        u_minus[i] = l[i];
        gamma = abs(s[i] + p);
        if(gamma < min_gamma){
          min_gamma = gamma;
          twist_index = i;
        }
        p -= shift;
      }
      else{
        if(is_nan(u_minus[i])){ //l[i]==0
          u_minus[i] = d[i];
        }
        //special case for next iteration of the for loop: p==inf, d_minus==inf, t==0
        i--;
        if(i < 0){
          break;
        }
        u_minus[i] = 0;
        gamma = abs(s[i] + d[i]);
        if(gamma < min_gamma){
          min_gamma = gamma;
          twist_index = i;
        }
        p = d[i] - shift;
      }
    }*/
    /*if(is_nan(p)){
      if(is_inf(d_minus)){ //since previous p==inf
        p = d[i] - shift;
        if(is_nan(t)){ //d[i]==0
          u_minus[i] = l[i];
        }
        else if(is_nan(u_minus[i])){ //t==inf && l[i]==0
          u_minus[i]=1;
        }
      }
    }*/
  }
  return twist_index;
}

/*
int getSturmCount(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiagSquared, double shift){
  int n = subdiagSquared.size();
  int count=0;
  double f1 = 1;
  double f2 = diag[0] - shift;
  count += f2>0;
  for(int i=0; i < n; i++){
    double f3 = (diag[i+1] - shift) * f2 - subdiagSquared[i] * f1;
    if(is_inf(f3)){
      return n+1;
    }
    count += (f3*f2>0) + (f2==0); //zero at the end does not count
    f1=f2;
    f2=f3;
  }
  return count;
}
 */
/*
int getSturmCountLdl(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiagSquared, double shift){
  double d = diag[0] - shift;
  int count = d>0;
  for(int i = 1; i < diag.size(); i++){
//    l = subdiag[i] / d_plus;
//    d_plus = diag[i+1] - shift - l * subdiag[i];
    d = diag[i] - shift - subdiagSquared[i-1] / d;
    count += d>0;
  }
  return count;
}
*/

int getSturmCountLdlLdl(const Eigen::VectorXd& d, const Eigen::VectorXd& l, double shift){
  int n = l.size();
  double s = -shift;
  double l_plus, d_plus;
  int count = 0;
  for(int i=0; i < n; i++){
    d_plus = s + d[i];
    count += d_plus>=0;
    //l_plus = l[i] * d[i] / d_plus;
    if (is_inf(d_plus) && is_inf(s)){ // this happens if d_plus==0 -> in next iteration d_plus==inf and s==inf
      s = l[i] * l[i] * d[i] - shift;
    }
    else{
      s = l[i] * l[i] * d[i] * s / d_plus - shift;
    }
  }
  d_plus = s + d[n];
  count += d_plus>=0;
  return count;
}

/*
void eigenvalsBisect(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag,  Eigen::VectorXd& eigenvals){
  int n = diag.size();
  double eps = 1e-11;

  double min_eigval = diag[0] - abs(subdiag[0]); //calc bounds on eigvals using Gresgorin discs
  double max_eigval = diag[0] + abs(subdiag[0]);
  for(int i=1;i<n-1;i++){
    min_eigval = std::min(min_eigval, diag[i] - abs(subdiag[i]) - abs(subdiag[i-1]));
    max_eigval = std::max(max_eigval, diag[i] + abs(subdiag[i]) + abs(subdiag[i-1]));
  }
  min_eigval = std::min(min_eigval, diag[n-1] - abs(subdiag[n-2]));
  max_eigval = std::max(max_eigval, diag[n-1] + abs(subdiag[n-2]));

  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();
  for(int i=0; i<n; i++){
    double low=min_eigval;
    double high=max_eigval;
    while(abs((high-low)/low)>eps){
      double mid=(high+low)*0.5;
      if(getSturmCountLdl(diag, subdiagSquared, mid)>i){
        low=mid;
      }
      else{
        high=mid;
      }
    }
    eigenvals[i]=(high+low)*0.5;
  }
}
 */

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

#define BISECT_K 8
/*
void eigenvalsBisect2(const Eigen::VectorXd& diag, const Eigen::VectorXd& subdiag,  Eigen::VectorXd& eigenvals){
  int n = diag.size();
  double eps = 1e-11;

  double min_eigval = diag[0] - abs(subdiag[0]); //calc bounds on eigvals using Gresgorin discs
  double max_eigval = diag[0] + abs(subdiag[0]);
  for(int i=1;i<n-1;i++){
    min_eigval = std::min(min_eigval, diag[i] - abs(subdiag[i]) - abs(subdiag[i-1]));
    max_eigval = std::max(max_eigval, diag[i] + abs(subdiag[i]) + abs(subdiag[i-1]));
  }
  min_eigval = std::min(min_eigval, diag[n-1] - abs(subdiag[n-2]));
  max_eigval = std::max(max_eigval, diag[n-1] + abs(subdiag[n-2]));

  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();
  Eigen::Array<double, BISECT_K, 1> low;
  Eigen::Array<double, BISECT_K, 1> high;
  for(int i=0; i<n; i += BISECT_K){
    low=Eigen::Array<double, BISECT_K, 1>::Constant(min_eigval);
    high=Eigen::Array<double, BISECT_K, 1>::Constant(max_eigval);
    while((abs((high-low)/low)>eps).any()){
      Eigen::Array<double, BISECT_K, 1> shifts = (high+low)*0.5;
      Eigen::Array<double, BISECT_K, 1> d = diag[0] - shifts;
      Eigen::Array<int, BISECT_K, 1> counts = (d>0).cast<int>();
      for(int j = 1; j < diag.size(); j++){ //Sturm count via LDL factorization
        d = diag[j] - shifts - subdiagSquared[j-1] / d;
        counts += (d>0).cast<int>();
      }
      for(int j=0;j<BISECT_K;j++){
        if(counts[j]>i+j){
          low[j]=shifts[j];
        }
        else{
          high[j]=shifts[j];
        }
      }
    }
    int n_valid = std::min(BISECT_K,n-i);
    eigenvals.segment(i,n_valid)=((high+low)*0.5).head(n_valid);
  }
}
 */


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
  double min_rel_sep=1e-1;
  double max_ele_growth = 1.7;
  int n=diag.size();
  Eigen::VectorXd high(n), low(n);
  double min_eigval;
  double max_eigval;
  getGresgorin(diag, subdiag, min_eigval, max_eigval);
  double span = max_eigval - min_eigval;
  Eigen::VectorXd subdiagSquared = subdiag.array() * subdiag.array();
  eigenvalsBisect3(diag, subdiagSquared, min_eigval, max_eigval, low, high);
  Eigen::VectorXd l(n-1), d(n);
  double shift0 = min_eigval;
  get_ldl(diag, subdiag, shift0, l, d);
  Eigen::VectorXd l2(n-1), d2(n), l3(n-1), d3(n), l_plus(n-1), d_plus(n), u_minus(n-1), d_minus(n), s(n), p(n);

  double T_norm = sqrt((diag.array()*diag.array()).sum() + subdiagSquared.sum());
  for(int i=0;i<n;i++){
    eigenvals[i]=(high[i]+low[i])*0.5;
    double shift = 0;
    int twist_idx;
    double low_gap = i==0 ? std::numeric_limits<double>::infinity() : low[i-1] - high[i];
    double high_gap = i==n-1 ? std::numeric_limits<double>::infinity() : low[i] - high[i+1];
    double min_gap = std::min(low_gap, high_gap);
    if(abs(low_gap / eigenvals[i]) > min_rel_sep || abs(high_gap / eigenvals[i]) < min_rel_sep){
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
      double min_element_growth = std::numeric_limits<double>::infinity();
      for(double sh : shifts){
        double element_growth = get_perturbed_shifted_ldl(d, l, sh, l3, d3);
        //cout << i << " element growth: " << element_growth << " at " << sh << endl;
        if(element_growth<min_element_growth){
          l2.swap(l3);
          d2.swap(d3);
          shift = sh + shift0;
          min_element_growth=element_growth;
        }
        if(element_growth<=max_ele_growth){
          break;
        }
      }

      //refine eigenvalue
      double shifted_low=low[i] - shift;
      double shifted_high=high[i] - shift;
      eigenvalBisectRefine(d2, l2, shifted_low, shifted_high, i);
      double shifted_eigenval = (shifted_low + shifted_high) * 0.5;
      twist_idx = get_twisted_factorization(d2, l2, shifted_eigenval, l_plus, u_minus);
    }
    else {
      twist_idx = get_twisted_factorization(d, l, eigenvals[i] - shift0, l_plus, u_minus);
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
    vec /= vec.norm();
  }
}

void symmetricEigenSolver(const Eigen::MatrixXd& A, Eigen::VectorXd& eigenvals, Eigen::MatrixXd& eigenvecs){
  Eigen::MatrixXd packed;
  block_householder_tridiag3(A, packed);
  Eigen::VectorXd diag = packed.diagonal();
  Eigen::VectorXd subdiag = packed.diagonal(1);
  mrrr(diag, subdiag, eigenvals, eigenvecs);
  //Eigen::MatrixXd q = Eigen::MatrixXd::Identity(eigenvecs.rows(), eigenvecs.cols());
  block_apply_packed_Q3(packed, eigenvecs);
  //eigenvecs = q * eigenvecs;
}


}
}

#endif //STAN_OPENCL
#endif //STAN_MATH_PRIM_MAT_FUN_OPENCL_EIGENDECOMPOSITION_HPP
