#ifndef SYMMETRC_EIGENSOLVER_CL_HPP
#define SYMMETRC_EIGENSOLVER_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/mrrr.hpp>
#include <stan/math/opencl/tridiagonalization.hpp>

stan::math::matrix_cl<double> A_glob_cl;

namespace stan {
namespace math {

matrix_cl<double> A2_glob_cl;
matrix_cl<double> packed_glob_cl;
matrix_cl<double> diag_glob;
matrix_cl<double> subdiag_glob;
matrix_cl<double> eigenvecs_glob;
matrix_cl<double> eigenvals_glob;
matrix_cl<double> eigenvecs2_glob;
matrix_cl<double> eigenvals2_glob;

template <bool need_eigenvectors = true>
void symmetric_eigensolver(const matrix_cl<double>& A,
                           matrix_cl<double>& eigenvalues,
                           matrix_cl<double>& eigenvectors) {
  int size = A.rows();
  matrix_cl<double> packed;
  if (is_first) {
      packed_glob_cl = packed;
      A2_glob_cl = A;
      std::cout
          << "A diff: "
          << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
          << std::endl;
  } else {
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;
    std::cout
        << "A2 diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A2_glob_cl))).maxCoeff()
        << std::endl;
  }
  internal::block_householder_tridiag_cl(A, packed);
  if (is_first) {
    packed_glob_cl = packed;
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;
  } else {

    std::cout << "packed diff: "
              << stan::math::from_matrix_cl(
                     max_2d(fabs(packed - packed_glob_cl)))
                     .maxCoeff()
              << std::endl;
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;
  }
  //  matrix_cl<double> packed2;
  //  internal::block_householder_tridiag_cl(A, packed2, 3);

  //  std::cout
  //          << "packed diff: "
  //          << stan::math::from_matrix_cl(max_2d(fabs(packed -
  //          packed2))).maxCoeff()
  //          << std::endl;
  matrix_cl<double> diag = diagonal(packed);
  matrix_cl<double> subdiag = diagonal(
      block_zero_based(packed, 0, 1, packed.rows() - 1, packed.cols() - 1));
//  matrix_cl<double> diag2 = diagonal(packed);
//  matrix_cl<double> subdiag2 = diagonal(
//      block_zero_based(packed2, 0, 1, packed.rows() - 1, packed.cols() - 1));

  if (is_first) {
    diag_glob = diag;
    subdiag_glob = subdiag;
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;
  } else {
    std::cout
        << "diag diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(diag - diag_glob))).maxCoeff()
        << std::endl;
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;

    std::cout << "subdiag diff: "
              << stan::math::from_matrix_cl(
                     max_2d(fabs(subdiag - subdiag_glob)))
                     .maxCoeff()
              << std::endl;
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;
  }
//  std::cout << "diag diff: "
//            << stan::math::from_matrix_cl(max_2d(fabs(diag - diag2))).maxCoeff()
//            << std::endl;
//  std::cout
//      << "subdiag diff: "
//      << stan::math::from_matrix_cl(max_2d(fabs(subdiag - subdiag2))).maxCoeff()
//      << std::endl;
  internal::tridiagonal_eigensolver_cl<need_eigenvectors>(
      diag, subdiag, eigenvalues, eigenvectors);
  if (is_first) {
    eigenvals_glob = eigenvalues;
    eigenvecs_glob = eigenvectors;
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;
  } else {
    std::cout
        << "eigenvalues diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(eigenvalues - eigenvals_glob))).maxCoeff()
        << std::endl;

    std::cout << "eigenvectors diff: "
              << stan::math::from_matrix_cl(
                     max_2d(fabs(eigenvectors - eigenvecs_glob)))
                     .maxCoeff()
              << std::endl;
    std::cout
        << "A diff: "
        << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
        << std::endl;
  }
//  std::cout << "diag diff after first: "
//            << stan::math::from_matrix_cl(max_2d(fabs(diag - diag2))).maxCoeff()
//            << std::endl;
//  std::cout
//      << "subdiag diff after first: "
//      << stan::math::from_matrix_cl(max_2d(fabs(subdiag - subdiag2))).maxCoeff()
//      << std::endl;
//  matrix_cl<double> eigenvalues2, eigenvectors2;
//  for (int i = 0; i < 20; i++) {
//    internal::tridiagonal_eigensolver_cl<need_eigenvectors>(
//        diag2, subdiag2, eigenvalues2, eigenvectors2);
//    std::cout << "eigenvalues diff: "
//              << stan::math::from_matrix_cl(
//                     max_2d(fabs(eigenvalues - eigenvalues2)))
//                     .maxCoeff()
//              << std::endl;
//    std::cout << "eigenvalues avg diff: "
//              << stan::math::from_matrix_cl(
//                     sum_2d(fabs(eigenvalues - eigenvalues2)))
//                         .sum()
//                     / A.rows()
//              << std::endl;
//  }
  if (need_eigenvectors) {
    //    std::cout << "eigenvectors diff: "
    //              << stan::math::from_matrix_cl(
    //                     max_2d(fabs(eigenvectors - eigenvectors2)))
    //                     .maxCoeff()
    //              << std::endl;
    //    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(A.rows(), A.cols());
    //    T.diagonal() = from_matrix_cl(diag);
    //    T.diagonal(-1) = T.diagonal(1) = from_matrix_cl(subdiag);
    //    std::cout << "vec ortho err: "
    //              << (from_matrix_cl(eigenvectors * transpose(eigenvectors))
    //                  - Eigen::MatrixXd::Identity(size, size))
    //                     .array()
    //                     .abs()
    //                     .maxCoeff()
    //              << std::endl;
    //    std::cout << "vec2 ortho err: "
    //              << (from_matrix_cl(eigenvectors2 * transpose(eigenvectors2))
    //                  - Eigen::MatrixXd::Identity(size, size))
    //                     .array()
    //                     .abs()
    //                     .maxCoeff()
    //              << std::endl;
    //    std::cout << "eigen eq err: "
    //              << ((T * from_matrix_cl(eigenvectors)
    //                   - from_matrix_cl(eigenvectors)
    //                         * from_matrix_cl(eigenvalues).asDiagonal())
    //                      .array()
    //                      .abs())
    //                     .array()
    //                     .abs()
    //                     .maxCoeff()
    //              << std::endl;
    //    std::cout << "eigen eq2 err: "
    //              << ((T * from_matrix_cl(eigenvectors2)
    //                   - from_matrix_cl(eigenvectors2)
    //                         * from_matrix_cl(eigenvalues2).asDiagonal())
    //                      .array()
    //                      .abs())
    //                     .array()
    //                     .abs()
    //                     .maxCoeff()
    //              << std::endl;
    //    std::cout << "eigen eq cross a err: "
    //              << ((T * from_matrix_cl(eigenvectors2)
    //                   - from_matrix_cl(eigenvectors2)
    //                         * from_matrix_cl(eigenvalues).asDiagonal())
    //                      .array()
    //                      .abs())
    //                     .array()
    //                     .abs()
    //                     .maxCoeff()
    //              << std::endl;
    //    std::cout << "eigen eq cross b err: "
    //              << ((T * from_matrix_cl(eigenvectors)
    //                   - from_matrix_cl(eigenvectors)
    //                         * from_matrix_cl(eigenvalues2).asDiagonal())
    //                      .array()
    //                      .abs())
    //                     .array()
    //                     .abs()
    //                     .maxCoeff()
    //              << std::endl;
    internal::block_apply_packed_Q_cl(packed, eigenvectors);
//    internal::block_apply_packed_Q_cl(packed2, eigenvectors2);

    if (is_first) {
      eigenvals2_glob = eigenvalues;
      eigenvecs2_glob = eigenvectors;
      is_first=false;
      std::cout
          << "A diff: "
          << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
          << std::endl;
    } else {
      std::cout
          << "eigenvalues2 diff: "
          << stan::math::from_matrix_cl(max_2d(fabs(eigenvalues - eigenvals2_glob))).maxCoeff()
          << std::endl;

      std::cout << "eigenvectors2 diff: "
                << stan::math::from_matrix_cl(
                       max_2d(fabs(eigenvectors - eigenvecs2_glob)))
                       .maxCoeff()
                << std::endl;
      std::cout
          << "A diff: "
          << stan::math::from_matrix_cl(max_2d(fabs(A - A_glob_cl))).maxCoeff()
          << std::endl;
    }
  }
}

}  // namespace math
}  // namespace stan

#endif
#endif
