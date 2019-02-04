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

//#define SKIP_CHECKS

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
#ifndef SKIP_CHECKS
  Mat R2 = R;
  R2.triangularView<Eigen::StrictlyLower>() = Mat::Constant(R.rows(), R.cols(), 0);
  cout << "R triang: " << (R - R2).array().abs().sum() << endl;
  cout << "reconstruct: " << (Q * R - a).array().abs().sum() << endl;
  cout << "ID: " << (Q.transpose() * Q).array().abs().sum() - Q.rows() << endl;
#endif
}

void chkEig(const Mat& a, const Mat& eigenvecs, const Vec& eigenvals){
#ifndef SKIP_CHECKS
  double A=a.rows();
  cout << "trace diff: " << a.diagonal().sum() - sum(eigenvals) << endl;
  cout << "eigen EQ: " << (a - eigenvecs * eigenvals.asDiagonal() * eigenvecs.transpose()).array().abs().sum() << endl;
  cout << "eigen total residuals: " << (a * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().abs().sum() << endl;
  cout << "eigen max residual: " << (a * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().abs().colwise().sum().maxCoeff() << endl;
  cout << "eigen total loss of orthogonality: " << (eigenvecs.transpose() * eigenvecs - Mat::Identity(A,A)).array().abs().sum() << endl;
  cout << "eigen max loss of orthogonality: " << (eigenvecs.transpose() * eigenvecs - Mat::Identity(A,A)).array().abs().maxCoeff() << endl;
#endif
}
/*
void chkTridiag(const Mat& a, Mat t, const Mat& q){
#ifndef SKIP_CHECKS
  cout << "reconstruct: " << (a - q * t * q.transpose()).array().abs().sum() << endl;
  //cout << q * t * q.transpose() << endl;
  cout << "ID: " << (q * q.transpose()).array().abs().sum() - q.rows() << endl;
  t.diagonal() = Mat::Constant(t.rows(), 1, 0);
  t.diagonal(1) = Mat::Constant(t.rows()-1, 1, 0);
  t.diagonal(-1) = Mat::Constant(t.rows()-1, 1, 0);
  cout << "tridiag: " << t.array().abs().sum() << endl;
#endif
}
 */

void chkTridiagPacked(const Mat& a, const Mat& packed){
#ifndef SKIP_CHECKS
  Mat t = Mat::Constant(a.rows(), a.cols(), 0);
  t.diagonal()=packed.diagonal();
  t.diagonal(1)=packed.diagonal(1);
  t.diagonal(-1)=packed.diagonal(1);
  if(a.rows()<10 && a.cols()<10){
    cout << "T: " << endl << t;
  }
  Mat q=Mat::Identity(a.rows(),a.cols());
  apply_packed_Q3(packed, q);
  if(a.rows()<10 && a.cols()<10){
    cout << "Q: " << endl << q;
  }
  cout << "reconstruct: " << (a - q * t * q.transpose()).array().abs().sum() << endl;
  cout << "ID: " << (q * q.transpose()).array().abs().sum() - q.rows() << endl;
  Mat a2=a;
  apply_packed_Q3(packed, a2);
  cout << "apply: " << (q * a - a2).array().abs().sum() << endl;
#endif
}
/*
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
  cout << "t" << endl << t << endl;
  cout << "q" << endl << q << endl;
  chkTridiag(a,t,q);
  householder_tridiag9(a, t, q);
  //block_householder_tridiag2(a, t, q, 2);
  cout << "t" << endl << t << endl;
  cout << "q" << endl << q << endl;
  chkTridiag(a,t,q);
}
 */
/*
int miniTest2() {
  Mat a(6, 6);
  a << 12, -51, 4, 4, 5, 9,
          6, 167, -68, 7, 9, 4,
          -6, 24, -41, 7, 6, 3,
          1, 2, 3, 4, 5, 4,
          1, 2, 3, 6, 5, 6,
          1, -7, 3, 4, 5, 6;
  a+=a.transpose().eval();
  Mat packed, vecs;
  Vec vals;

  cout << "a:" << endl << a << endl;
  householder_tridiag_packed(a, packed);
  cout << "packed:" << endl << packed << endl;
  chkTridiagPacked(a,packed);
}
*/

int main() {
  auto start = std::chrono::steady_clock::now();
  /*miniTest();
  miniTest2();
  return 0;*/
  cout.precision(4);
  int A = 2000;
  for(unsigned int i=443;i<1e7;i++) {
    cout << "i=" << i << endl;
    srand(i);

    Vec diag = Vec::Random(A).array();
    Vec subdiag = Vec::Random(A - 1).array();
    /*subdiag[2]=0;
    subdiag[3]=0;
    diag[5]=0;
    diag[6]=0;
    diag[10]=0;
    diag[11]=0;
    subdiag[10]=0;*/
    Vec eigenvals(A), d(A), d_plus(A), d_minus(A), l(A - 1), l_plus(A - 1), u(A - 1), u_minus(A - 1), s(A), p_(A), low(A), high(A);
    Mat eigenvecs(A, A), eigenvecs2(A, A);
    Mat T = Mat::Constant(A, A, 0);
    T.diagonal() = diag;
    T.diagonal(1) = subdiag;
    T.diagonal(-1) = subdiag;

    start = std::chrono::steady_clock::now();
    mrrr(diag, subdiag, eigenvals, eigenvecs2);
    if(A>=15) {
      cout << "mrrr: "
           << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
           << "ms" << endl;
    }
    eigenvecs = eigenvecs2.array().rowwise() / eigenvecs2.colwise().norm().array();
    if(eigenvecs2.unaryExpr([](double x){return is_nan(x);}).any()){
      cout << "nan in eigenvecs!" << endl;
    }
    if(eigenvecs2.unaryExpr([](double x){return is_inf(x);}).any()){
      cout << "inf in eigenvecs!" << endl;
    }
    /*if((eigenvecs.transpose() * eigenvecs - Mat::Identity(A,A)).array().maxCoeff()<0.1 ||
    (T * eigenvecs - eigenvecs * eigenvals.asDiagonal()).array().colwise().sum().maxCoeff() <0.1){
      continue;
    }*/
//    cout << endl;
    if (A < 15) {
      cout << "mat:" << endl;
      cout << diag << endl << endl;
      cout << subdiag << endl << endl;
      cout << "my results:" << endl;
      cout << "eigenvectors:" << endl << eigenvecs2 << endl << endl;
      cout << "eigenvectors:" << endl << eigenvecs << endl << endl;
      cout << "eigenvalues:" << endl << eigenvals << endl << endl;
//    cout << (T*eigenvecs.col(0)).eval() << endl << endl;
//    cout << (eigenvals[0]*eigenvecs.col(0)).eval() << endl << endl;
//      cout << T * eigenvecs << endl << endl;
//      cout << eigenvecs * eigenvals.asDiagonal() << endl << endl;
//      cout << eigenvecs.transpose() * eigenvecs;
    }
    /**/
    chkEig(T, eigenvecs, eigenvals);

    SelfAdjointEigenSolver<Mat> slv(T);
    cout << "CPU: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
         << "ms" << endl;
    eigenvals=slv.eigenvalues();
    eigenvecs=slv.eigenvectors();
    cout << "CPU total: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
         << "ms" << endl;
    if (A < 15) {
      cout << "eigen results:" << endl;
      cout << "eigenvectors:" << endl << eigenvecs.rowwise().reverse() << endl << endl;
      cout << "eigenvalues:" << endl << eigenvals.reverse() << endl << endl;
    }
    chkEig(T,eigenvecs,eigenvals);
    //break;
  }


//  Mat L=Mat::Identity(A,A), U=Mat::Identity(A,A);
//  L.diagonal(-1)=l;
//  U.diagonal(1)=u;
//  start = std::chrono::steady_clock::now();
//  get_ldl(diag,subdiag,l,d);
//  get_udu(diag,subdiag,u,d_minus);
//  cout << "twisted decomp: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  cout << "LDL: " << (T-L*d.asDiagonal()*L.transpose()).array().abs().sum() << endl;
//  cout << "UDU: " << (T-U*d_minus.asDiagonal()*U.transpose() ).array().abs().sum() << endl;


//  double shift=2;
//  get_shifted_ldl(d,l, shift, l_plus, d_plus, s);
//  get_shifted_udu(d,l, shift, u_minus, d_minus, p);
//  L.diagonal(-1)=l_plus;
//  U.diagonal(1)=u_minus;
//  cout << "shifted LDL: " << (T-(Mat::Identity(A,A).array()*shift).matrix()-L*d_plus.asDiagonal() *L.transpose()).array().abs().sum() << endl;
//  cout << "shifted UDU: " << (T-(Mat::Identity(A,A).array()*shift).matrix()-U*d_minus.asDiagonal()*U.transpose()).array().abs().sum() << endl;
//
//  get_perturbed_shifted_ldl(d,l, shift, l_plus, d_plus);
//  L.diagonal(-1)=l_plus;
//  cout << "perturbed shifted LDL: " << (T-(Mat::Identity(A,A).array()*shift).matrix()-L*d_plus.asDiagonal() *L.transpose()).array().abs().sum() << endl;



//  start = std::chrono::steady_clock::now();
//  dqds(diag, subdiag, eigenvals);
//  cout << "dqds: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  if (A < 10) {
//    cout << "eigvals: " << eigenvals << endl;
//  }
//  cout << "trace diff: " << sum(diag)-sum(eigenvals) << endl;
//
//
//  Vec subdiagSq=subdiag.array()*subdiag.array();
//  double min_eigval;
//  double max_eigval;
//  getGresgorin(diag, subdiag, min_eigval, max_eigval);
//  start = std::chrono::steady_clock::now();
//  eigenvalsBisect3(diag, subdiagSq, min_eigval, max_eigval, low, high);
//  cout << "bisect: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  eigenvals=(low+high)/2;
//  //cout << "eigvals: " << eigenvals << endl;
//  cout << "trace diff: " << sum(diag)-sum(eigenvals) << endl;


//  for(double shift  : {-1.2,-0.8,-0.300001,-0.3, -0.299999,-0.0001,0.0001,0.1,0.4,1.6,1.8}){
//    int count=getSturmCountLdl(diag,subdiagSq,shift);
//    int count2=0;
//    for(int i=0;i<eigenvals.size();i++){
//      count2 += eigenvals[i]>shift;
//    }
//    cout << "shift: " << shift << " Sturm: " << count << " actual: " << count2 << endl;
//  }

//  subdiag.array()-=500;
//  start = std::chrono::steady_clock::now();
//  dqds(diag.array(), subdiag, eigenvals);
//  cout << "dqds-500: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  //cout << "eigvals: " << eigenvals << endl;
//  cout << "trace diff: " << sum(diag)-sum(eigenvals) << endl;
//
  return 0;

/*
  Mat a = Mat::Random(A, A);
  a+=a.transpose().eval();
  a=-100*a;
  Mat t,q, vecs, packed, qa, qa2, q2;
  Vec vals;

//  SelfAdjointEigenSolver<Mat> slv(a);
//  cout << "CPU: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  vals=slv.eigenvalues();
//  vecs=slv.eigenvectors();
//  cout << "CPU total: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkEig(a,vecs,vals);


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

//  start = std::chrono::steady_clock::now();
//  householder_tridiag(a, t, q);
//  cout << "\t\tCPU my basic: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//  start = std::chrono::steady_clock::now();
//  householder_tridiag2(a, t, q);
//  cout << "\t\tCPU my basic2: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//  start = std::chrono::steady_clock::now();
//  householder_tridiag3(a, t, q);
//  cout << "\t\tCPU my basic3: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);



//  start = std::chrono::steady_clock::now();
//  householder_tridiag4(a, t, q);
//  cout << "\t\tCPU my basic4: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//
//  start = std::chrono::steady_clock::now();
//  householder_tridiag5(a, t, q);
//  cout << "\t\tCPU my basic5: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//
//  start = std::chrono::steady_clock::now();
//  householder_tridiag6(a, t, q);
//  cout << "\t\tCPU my basic6: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//
//  start = std::chrono::steady_clock::now();
//  householder_tridiag7(a, t, q);
//  cout << "\t\tCPU my basic7: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//
//  start = std::chrono::steady_clock::now();
//  householder_tridiag8(a, t, q);
//  cout << "\t\tCPU my basic8: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//
//  start = std::chrono::steady_clock::now();
//  householder_tridiag9(a, t, q);
//  cout << "\t\tCPU my basic9: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);


  start = std::chrono::steady_clock::now();
  householder_tridiag_packed(a, packed);
  cout << "\t\tCPU my basic packed: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;

  start = std::chrono::steady_clock::now();
  block_householder_tridiag3(a, packed);
  cout << "\t\tCPU my blocked packed: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;

  t = Mat::Constant(a.rows(), a.cols(), 0);
  t.diagonal()=packed.diagonal();
  t.diagonal(1)=packed.diagonal(1);
  t.diagonal(-1)=packed.diagonal(1);

//  start = std::chrono::steady_clock::now();
//  q=Mat::Identity(a.rows(),a.cols());
//  apply_packed_Q(packed,q);
//  cout << "\t\tCPU apply packed Q: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//
//  start = std::chrono::steady_clock::now();
//  q=Mat::Identity(a.rows(),a.cols());
//  apply_packed_Q2(packed,q);
//  cout << "\t\tCPU apply packed Q2: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);

  start = std::chrono::steady_clock::now();
  q=Mat::Identity(a.rows(),a.cols());
  apply_packed_Q3(packed,q);
  cout << "\t\tCPU apply packed Q3: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
  qa2=a;
  apply_packed_Q3(packed,qa2);
  q2=q;


  start = std::chrono::steady_clock::now();
  q=Mat::Identity(a.rows(),a.cols());
  apply_packed_Q4(packed,q);
  cout << "\t\tCPU apply packed Q4: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
  qa=a;
  apply_packed_Q4(packed,qa);
  cout << "qa: " << (qa-qa2).array().abs().sum() << endl;
  cout << "q: " << (q-q2).array().abs().sum() << endl;


  start = std::chrono::steady_clock::now();
  q=Mat::Identity(a.rows(),a.cols());
  block_apply_packed_Q(packed,q);
  cout << "\t\tCPU block apply packed Q: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
  qa=a;
  block_apply_packed_Q(packed,qa);
  cout << "qa: " << (qa-qa2).array().abs().sum() << endl;
  cout << "q: " << (q-q2).array().abs().sum() << endl;

  start = std::chrono::steady_clock::now();
  q=Mat::Identity(a.rows(),a.cols());
  block_apply_packed_Q2(packed,q);
  cout << "\t\tCPU block apply packed Q2: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  chkTridiag(a,t,q);
  qa=a;
  block_apply_packed_Q2(packed,qa);
  cout << "qa: " << (qa-qa2).array().abs().sum() << endl;
  cout << "q: " << (q-q2).array().abs().sum() << endl;

//  start = std::chrono::steady_clock::now();
//  q=Mat::Identity(a.rows(),a.cols());
//  block_apply_packed_Q3(packed,q,3);
//  cout << "\t\tCPU block apply packed Q3: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chkTridiag(a,t,q);
//  qa=a;
//  block_apply_packed_Q3(packed,qa,3);
//  cout << "qa: " << (qa-qa2).array().abs().sum() << endl;

  start = std::chrono::steady_clock::now();
  qa=q*a;
  cout << "\t\tCPU multiply q*a: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  cout << "qa: " << (qa-qa2).array().abs().sum() << endl;

  start = std::chrono::steady_clock::now();
  qa=a*q;
  cout << "\t\tCPU multiply a*q: "
       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
       << "ms" << endl;
  cout << "aq: " << (qa-qa2).array().abs().sum() << endl;
//  cout << q.transpose()*a*q << endl;
//  cout << q*a*q.transpose() << endl;
//  cout << endl << t << endl;


//  for(int b=5;b<400;b+=5) {
//    cout << "b = " << b << endl;
//
////    cout << "b = " << b << endl;
////    start = std::chrono::steady_clock::now();
////    block_householder_tridiag2(a, t, q, b);
////    cout << "\t\tCPU my blocked: "
////         << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
////         << "ms" << endl;
////    chkTridiag(a, t, q);
//  }


//  //force kernel compilation
//  cl::Kernel kernel_1 = opencl_context.get_kernel("householder_QR_1");
//  cl::Kernel kernel_2 = opencl_context.get_kernel("matrix_multiply");
//
//  start = std::chrono::steady_clock::now();
//  block_householder_qr_gpu_hybrid(a, Q, R, 120);
//  cout << "GPU - hybrid my block 120: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chk(a,Q,R);
//
//  block_householder_qr_gpu(a, Q, R, 160);
//  start = std::chrono::steady_clock::now();
//  block_householder_qr_gpu(a, Q, R, 160);
//  cout << "GPU my block 160: "
//       << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()
//       << "ms" << endl;
//  chk(a,Q,R);

  return 0;*/
}