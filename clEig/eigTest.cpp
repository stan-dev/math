//
// Created by tadej on 7. 12. 2018.
//

#include <iostream>

#define STAN_OPENCL
#define OPENCL_PLATFORM_ID 0
#define OPENCL_DEVICE_ID 0

#include <stan/math.hpp>
#include <stan/math/gpu/add.hpp>
#include <stan/math/gpu/copy.hpp>
#include <chrono>
#include <Eigen/QR>
#include <stan/math/gpu/eigendecomposition.hpp>


using namespace Eigen;
using namespace std;
using namespace stan::math;

using Mat=Matrix<double, Dynamic, Dynamic>;
using Vec=VectorXd;

void p(const Eigen::MatrixXd& a) {
  std::cout << a << std::endl;
}

void s(const Eigen::MatrixXd& a) {
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

void chk(const Mat& a, const Mat& Q, const Mat& R){
  Mat R2 = R;
  R2.triangularView<Eigen::StrictlyLower>() = Mat::Constant(R.rows(), R.cols(), 0);
  cout << "R triang: " << (R - R2).array().abs().sum() << endl;
  cout << "reconstruct: " << (Q * R - a).array().abs().sum() << endl;
  cout << "ID: " << (Q.transpose() * Q).array().abs().sum() - Q.rows() << endl;
}

void chkEig(const Mat& a, const Mat& vecs, const Vec& vals){
  cout << "eig: " << (a*vecs - vecs*vals.asDiagonal()).array().abs().sum() << endl;
}

void chkTridiag(const Mat& a, Mat t, const Mat& q){
  cout << "reconstruct: " << (a - q * t * q.transpose()).array().abs().sum() << endl;
  //cout << q * t * q.transpose() << endl;
  cout << "ID: " << (q * q.transpose()).array().abs().sum() - q.rows() << endl;
  t.diagonal() = Mat::Constant(t.rows(), 1, 0);
  t.diagonal(1) = Mat::Constant(t.rows()-1, 1, 0);
  t.diagonal(-1) = Mat::Constant(t.rows()-1, 1, 0);
  cout << "tridiag: " << t.array().abs().sum() << endl;
}

int miniTest() {
  Mat a(5, 5);
  a << 12, -51, 4, 4, 5,
        6, 167, -68, 7, 9,
        -6, 24, -41, 7, 6,
        1, 2, 3, 4, 5,
        1, 2, 3, 4, 5;
  a+=a.transpose().eval();
  Mat t,q, vecs;
  Vec vals;

  cout << "a:" << endl << a << endl;
  householder_tridiag2(a, t, q);
  //block_householder_tridiag(a, t, q, 2);
  chkTridiag(a,t,q);
}

int main() {
  //miniTest();
  //return 0;

  int A = 601;
  const int MAX_BLOCK=400;
  Mat a = Mat::Random(A, A);
  a+=a.transpose().eval();
  Mat t,q, vecs;
  Vec vals;

  auto start = std::chrono::steady_clock::now();
  /*SelfAdjointEigenSolver<Mat> slv(a);
  cout << "CPU: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  vals=slv.eigenvalues();
  vecs=slv.eigenvectors();
  cout << "CPU total: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkEig(a,vecs,vals);
*/

  start = std::chrono::steady_clock::now();
  Tridiagonalization<Mat> slv2(a);
  cout << "CPU: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  t=slv2.matrixT();
  q=slv2.matrixQ();
  cout << "CPU total: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);

  start = std::chrono::steady_clock::now();
  householder_tridiag(a, t, q);
  cout << "CPU my basic: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);

  start = std::chrono::steady_clock::now();
  householder_tridiag2(a, t, q);
  cout << "CPU my basic2: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);

/*
  start = std::chrono::steady_clock::now();
  block_householder_tridiag(a, t, q);
  cout << "CPU my blocked: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
*/
  /*
  //force kernel compilation
  cl::Kernel kernel_1 = opencl_context.get_kernel("householder_QR_1");
  cl::Kernel kernel_2 = opencl_context.get_kernel("matrix_multiply");

  start = std::chrono::steady_clock::now();
  block_householder_qr_gpu_hybrid(a, Q, R, 120);
  cout << "GPU - hybrid my block 120: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chk(a,Q,R);

  block_householder_qr_gpu(a, Q, R, 160);
  start = std::chrono::steady_clock::now();
  block_householder_qr_gpu(a, Q, R, 160);
  cout << "GPU my block 160: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chk(a,Q,R);
*/
  return 0;
}