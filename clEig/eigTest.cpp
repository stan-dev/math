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
  Mat a(6, 6);
  a << 12, -51, 4, 4, 5, 9,
        6, 167, -68, 7, 9, 4,
        -6, 24, -41, 7, 6, 3,
        1, 2, 3, 4, 5, 4,
        1, 2, 3, 6, 5, 6,
        1, -7, 3, 4, 5, 6;
  a+=a.transpose().eval();
  Mat t,q, vecs;
  Vec vals;

  cout << "a:" << endl << a << endl;
  householder_tridiag8(a, t, q);
  cout << "t" << endl;
  cout << t << endl;
  cout << "q" << endl;
  cout << q << endl;
  chkTridiag(a,t,q);
  householder_tridiag9(a, t, q);
  //block_householder_tridiag2(a, t, q, 2);
  cout << "t" << endl;
  cout << t << endl;
  cout << "q" << endl;
  cout << q << endl;
  chkTridiag(a,t,q);
}

int main() {
  //miniTest();
  //return 0;

  int A = 1000;
  const int MAX_BLOCK=170;
  Mat a = Mat::Random(A, A);
  a+=a.transpose().eval();
  a=-a;
  Mat t,q, vecs;
  Vec vals;
  auto start = std::chrono::steady_clock::now();
/*
  SelfAdjointEigenSolver<Mat> slv(a);
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
  cout << "\t\tCPU total: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
/*
  start = std::chrono::steady_clock::now();
  householder_tridiag(a, t, q);
  cout << "\t\tCPU my basic: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);

  start = std::chrono::steady_clock::now();
  householder_tridiag2(a, t, q);
  cout << "\t\tCPU my basic2: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);

  start = std::chrono::steady_clock::now();
  householder_tridiag3(a, t, q);
  cout << "\t\tCPU my basic3: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
  */

/*
  start = std::chrono::steady_clock::now();
  householder_tridiag4(a, t, q);
  cout << "\t\tCPU my basic4: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
  */


  start = std::chrono::steady_clock::now();
  householder_tridiag5(a, t, q);
  cout << "\t\tCPU my basic5: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);


  start = std::chrono::steady_clock::now();
  householder_tridiag6(a, t, q);
  cout << "\t\tCPU my basic6: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);


  start = std::chrono::steady_clock::now();
  householder_tridiag7(a, t, q);
  cout << "\t\tCPU my basic7: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);


  start = std::chrono::steady_clock::now();
  householder_tridiag8(a, t, q);
  cout << "\t\tCPU my basic8: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);


  start = std::chrono::steady_clock::now();
  householder_tridiag9(a, t, q);
  cout << "\t\tCPU my basic9: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);

/*
  for(int b=10;b<35;b+=1) {
    cout << "b = " << b << endl;
    start = std::chrono::steady_clock::now();
    block_householder_tridiag2(a, t, q, b);
    cout << "\t\tCPU my blocked: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
         << "ms" << endl;
    chkTridiag(a, t, q);
  }
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