#ifndef STAN_MATH_PRIM_FUN_MATRIX_EXP_ACTION_HANDLER_HPP
#define STAN_MATH_PRIM_FUN_MATRIX_EXP_ACTION_HANDLER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/ceil.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <cmath>

namespace stan {
namespace math {

/**
 * The implementation of the work by Awad H. Al-Mohy and Nicholas J. Higham
 * "Computing the Action of the Matrix Exponential,
 * with an Application to Exponential Integrators"
 * Read More: https://epubs.siam.org/doi/abs/10.1137/100788860
 * See also:
 * https://www.mathworks.com/matlabcentral/fileexchange/29576-matrix-exponential-times-a-vector
 *
 * Calculates exp(mat*t)*b, where mat & b are matrices,
 * and t is double.
 */
class matrix_exp_action_handler {
  const int _p_max = 8;
  const int _m_max = 55;
  const double _tol = 1.1e-16;  // from the paper, double precision: 2^-53
  const std::vector<double> _theta_m{
      2.22044605e-16, 2.58095680e-08, 1.38634787e-05, 3.39716884e-04,
      2.40087636e-03, 9.06565641e-03, 2.38445553e-02, 4.99122887e-02,
      8.95776020e-02, 1.44182976e-01, 2.14235807e-01, 2.99615891e-01,
      3.99777534e-01, 5.13914694e-01, 6.41083523e-01, 7.80287426e-01,
      9.30532846e-01, 1.09086372e+00, 1.26038106e+00, 1.43825260e+00,
      1.62371595e+00, 1.81607782e+00, 2.01471078e+00, 2.21904887e+00,
      2.42858252e+00, 2.64285346e+00, 2.86144963e+00, 3.08400054e+00,
      3.31017284e+00, 3.53966635e+00, 3.77221050e+00, 4.00756109e+00,
      4.24549744e+00, 4.48581986e+00, 4.72834735e+00, 4.97291563e+00,
      5.21937537e+00, 5.46759063e+00, 5.71743745e+00, 5.96880263e+00,
      6.22158266e+00, 6.47568274e+00, 6.73101590e+00, 6.98750228e+00,
      7.24506843e+00, 7.50364669e+00, 7.76317466e+00, 8.02359473e+00,
      8.28485363e+00, 8.54690205e+00, 8.80969427e+00, 9.07318789e+00,
      9.33734351e+00, 9.60212447e+00, 9.86749668e+00, 1.01334283e+01,
      1.03998897e+01, 1.06668532e+01, 1.09342929e+01, 1.12021845e+01,
      1.14705053e+01, 1.17392341e+01, 1.20083509e+01, 1.22778370e+01,
      1.25476748e+01, 1.28178476e+01, 1.30883399e+01, 1.33591369e+01,
      1.36302250e+01, 1.39015909e+01, 1.41732223e+01, 1.44451076e+01,
      1.47172357e+01, 1.49895963e+01, 1.52621795e+01, 1.55349758e+01,
      1.58079765e+01, 1.60811732e+01, 1.63545578e+01, 1.66281227e+01,
      1.69018609e+01, 1.71757655e+01, 1.74498298e+01, 1.77240478e+01,
      1.79984136e+01, 1.82729215e+01, 1.85475662e+01, 1.88223426e+01,
      1.90972458e+01, 1.93722711e+01, 1.96474142e+01, 1.99226707e+01,
      2.01980367e+01, 2.04735082e+01, 2.07490816e+01, 2.10247533e+01,
      2.13005199e+01, 2.15763782e+01, 2.18523250e+01, 2.21283574e+01};

 public:
  /**
   * Perform the matrix exponential action exp(A*t)*B
   * @param [in] mat matrix A
   * @param [in] b matrix B
   * @param [in] t double t, e.g. time.
   * @return matrix exp(A*t)*B
   */
  template <typename EigMat1, typename EigMat2,
            require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
            require_all_st_same<double, EigMat1, EigMat2>* = nullptr>
  inline Eigen::MatrixXd action(const EigMat1& mat, const EigMat2& b,
                                const double& t = 1.0) {
    Eigen::MatrixXd A = mat;
    double mu = A.trace() / A.rows();
    A.diagonal().array() -= mu;

    int m, s;
    set_approx_order(A, b, t, m, s);

    double eta = exp(t * mu / s);

    const auto& b_eval = b.eval();
    Eigen::MatrixXd f = b_eval;
    Eigen::MatrixXd bi = b_eval;

    for (int i = 0; i < s; ++i) {
      //  maximum absolute row sum, aka inf-norm of matrix operator
      double c1 = matrix_operator_inf_norm(bi);
      for (int j = 0; j < m; ++j) {
        bi = (t / (s * (j + 1))) * (A * bi);
        f += bi;
        double c2 = matrix_operator_inf_norm(bi);
        if (c1 + c2 < _tol * matrix_operator_inf_norm(f))
          break;
        c1 = c2;
      }
      f *= eta;
      bi = f;
    }

    return f;
  }

  /**
   * Eigen expression for matrix operator infinity norm
   *
   * @param x matrix
   */
  double matrix_operator_inf_norm(Eigen::MatrixXd const& x) {
    return x.cwiseAbs().rowwise().sum().maxCoeff();
  }

  /**
   * Estimate the 1-norm of mat^m
   *
   * See A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
   *   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
   *   970-989, 2009.
   *
   * For positive matrices the results is exact. Otherwise it falls
   * back to Eigen's norm, which is only
   * efficient for small & medium-size matrices (n < 100). Large size
   * matrices require a more efficient 1-norm approximation
   * algorithm such as normest1. See, e.g.,
   * https://hg.savannah.gnu.org/hgweb/octave/file/e35866e6a2e0/scripts/linear-algebra/normest1.m
   *
   * @param mat matrix
   * @param m power
   */
  template <typename EigMat1, require_all_eigen_t<EigMat1>* = nullptr,
            require_all_st_same<double, EigMat1>* = nullptr>
  double mat_power_1_norm(const EigMat1& mat, int m) {
    if ((mat.array() > 0.0).all()) {
      Eigen::VectorXd e = Eigen::VectorXd::Constant(mat.rows(), 1.0);
      for (int j = 0; j < m; ++j) {
        e = mat.transpose() * e;
      }
      return e.lpNorm<Eigen::Infinity>();
    } else {
      return mat.pow(m).cwiseAbs().colwise().sum().maxCoeff();
    }
  }

  /**
   * Approximation is based on parameter "m" and "s",
   * proposed in CODE FRAGMENT 3.1 of the reference. The
   * parameters depend on matrix A, as well as limits
   * "m_max" & "p_max", defined in table 3.1 and eq. 3.11,
   * respectively. Ideally one can choose "m_max" and
   * "p_max" to suite the specific computing precision needs,
   * here we just use the one suggested in the paper
   * (paragraph between eq. 3.11 & eq. 3.12).
   *
   * @param [in] mat matrix A
   * @param [in] b matrix B
   * @param [in] t double t in exp(A*t)*B
   * @param [out] m int parameter m; m >= 0, m < _m_max
   * @param [out] s int parameter s; s >= 1
   */
  template <typename EigMat1, typename EigMat2,
            require_all_eigen_t<EigMat1, EigMat2>* = nullptr,
            require_all_st_same<double, EigMat1, EigMat2>* = nullptr>
  inline void set_approx_order(const EigMat1& mat, const EigMat2& b,
                               const double& t, int& m, int& s) {
    if (t < _tol) {
      m = 0;
      s = 1;
      return;
    }

    // L1 norm
    double normA = mat.colwise().template lpNorm<1>().maxCoeff();
    if (normA < _tol) {
      m = 0;
      s = 1;
      return;
    }

    Eigen::VectorXd alpha(_p_max - 1);
    if (normA * (_m_max * b.cols())
        < 4.0 * _theta_m[_m_max] * _p_max * (_p_max + 3)) {
      alpha = Eigen::VectorXd::Constant(_p_max - 1, normA);
    } else {
      Eigen::VectorXd eta(_p_max);
      for (int p = 0; p < _p_max; ++p) {
        eta[p] = std::pow(mat_power_1_norm(mat, p + 2), 1.0 / (p + 2));
      }
      for (int p = 0; p < _p_max - 1; ++p) {
        alpha[p] = std::max(eta[p], eta[p + 1]);
      }
    }

    Eigen::MatrixXd mt = Eigen::MatrixXd::Zero(_p_max - 1, _m_max);
    for (int p = 1; p < _p_max; ++p) {
      for (int i = p * (p + 1) - 2; i < _m_max; ++i) {
        mt(p - 1, i) = alpha[p - 1] / _theta_m[i];
      }
    }

    Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(_m_max, 1, _m_max);
    Eigen::MatrixXd c = stan::math::ceil(mt) * u.asDiagonal();
    for (Eigen::MatrixXd::Index i = 0; i < c.size(); ++i) {
      if (c(i) <= 1.e-16) {
        c(i) = std::numeric_limits<double>::infinity();
      }
    }
    int cost = c.colwise().minCoeff().minCoeff(&m);
    if (std::isinf(cost)) {
      s = 1;
      return;
    }
    s = std::max(cost / m, 1);
  }
};

}  // namespace math
}  // namespace stan

#endif
