#ifndef STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_GLM_LPMF_HPP
#define STAN_MATH_PRIM_MAT_PROB_ORDERED_LOGISTIC_GLM_LPMF_HPP

#include <stan/math/prim/scal/meta/is_constant_struct.hpp>
#include <stan/math/prim/scal/meta/partials_return_type.hpp>
#include <stan/math/prim/scal/meta/operands_and_partials.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/arr/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/meta/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/scal/meta/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/mat/err/check_ordered.hpp>
#include <stan/math/prim/scal/fun/log1p_exp.hpp>
#include <cmath>

namespace stan {
namespace math {

void s(const Eigen::MatrixXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}
void s(const Eigen::VectorXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}
void s(const Eigen::RowVectorXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}

void s(const Eigen::ArrayXXd& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}
void s(const Eigen::Array<double,Eigen::Dynamic,1>& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
}
void s(const Eigen::Array<double,1,Eigen::Dynamic>& a) {
  std::cout << "(" << a.rows() << ", " << a.cols() << std::endl;
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

void p(const Eigen::ArrayXXd& a) {
  s(a);
  std::cout << a << std::endl;
}
void p(const Eigen::Array<double,Eigen::Dynamic,1>& a) {
  s(a);
  std::cout << a << std::endl;
}
void p(const Eigen::Array<double,1,Eigen::Dynamic>& a) {
  s(a);
  std::cout << a << std::endl;
}

#ifdef STAN_OPENCL
void p(const matrix_cl& a) {
  Eigen::MatrixXd b(a.rows(), a.cols());
  b = from_matrix_cl(a);
  std::cout << b << std::endl;
}
#endif

template<bool propto, typename T_y, typename T_x, typename T_beta, typename T_cuts>
typename stan::return_type<T_x, T_beta, T_cuts>::type
ordered_logistic_glm_lpmf(
        const T_y& y, const Eigen::Matrix<T_x, Eigen::Dynamic, Eigen::Dynamic>& x,
        const T_beta& beta, const T_cuts& cuts) {
  static const char* function = "ordered_logistic_glm_lpmf";
  typedef typename partials_return_type<T_y, T_x, T_beta, T_cuts>::type
          T_partials_return;

  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::Array;
  using Eigen::VectorXd;
  using std::exp;

  const size_t N_instances = x.rows();
  const size_t N_attributes = x.cols();
  const size_t N_classes = length(cuts) + 1;

  check_consistent_size(function, "Vector of dependent variables", y, N_instances);
  check_consistent_size(function, "Weight vector", beta, N_attributes);
  check_bounded(function, "Vector of dependent variables", y, 1, N_classes);
  check_ordered(function, "Cut-points", cuts);
  check_finite(function, "Final cut-point", cuts[N_classes - 2]);
  check_finite(function, "First cut-point", cuts[0]);

  if (size_zero(y, x, beta))
    return 0;

  if (!include_summand<propto, T_x, T_beta, T_cuts>::value)
    return 0;

  T_partials_return logp(0);
  const auto& x_val = value_of_rec(x);
  const auto& beta_val = value_of_rec(beta);
  const auto& cuts_val = value_of_rec(cuts);

  const auto& y_vec = as_column_vector_or_scalar(y);
  const auto& beta_val_vec = as_column_vector_or_scalar(beta_val);
  const auto& cuts_val_vec = as_column_vector_or_scalar(cuts_val);

  Array<double, Dynamic, 1> cuts_y1(N_instances), cuts_y2(N_instances);
  for (int i = 0; i < N_instances; i++) {
    if (y[i] != N_classes) {
      cuts_y1[i] = cuts_val_vec[y[i] - 1]; //y==C: INF
    }
    else{
      cuts_y1[i] = INFINITY;
    }
    if (y[i] != 1) {
      cuts_y2[i] = cuts_val_vec[y[i] - 2]; //y==1: -INF
    }
    else{
      cuts_y2[i] = -INFINITY;
    }
  }

  Array<double, Dynamic, 1> location = x_val * beta_val_vec;
  if(!location.allFinite()){
    check_finite(function, "Weight vector", beta);
    check_finite(function, "Matrix of independent variables", x);
  }

  Array<double, Dynamic, 1> cut2 = location - cuts_y2; //y==1: INF
  Array<double, Dynamic, 1> cut1 = location - cuts_y1; //y==C: -INF

  //TODO eliminate -
  //TODO select?
  Array<double, Dynamic, 1> m_log_1p_exp_cut1 = -cut1 * (cut1 > 0.0).cast<double>() - (-cut1.abs()).exp().log1p();
  Array<double, Dynamic, 1> m_log_1p_exp_m_cut2 = cut2 * (cut2 <= 0.0).cast<double>() - (-cut2.abs()).exp().log1p();
  //TODO common subexpressions or move selects up?
  logp = y.cwiseEqual(1).select(m_log_1p_exp_cut1,
                                y.cwiseEqual(N_classes).select(m_log_1p_exp_m_cut2,
                                                       m_log_1p_exp_m_cut2 + log1m_exp(cut1 - cut2).array() + m_log_1p_exp_cut1
                                )).sum();

  operands_and_partials<Matrix<T_x, Dynamic, Dynamic>, T_beta, T_cuts> ops_partials(x, beta, cuts);
  if (!is_constant_struct<T_x>::value && !is_constant_struct<T_beta>::value && !is_constant_struct<T_cuts>::value) {
    //TODO single fraction?
    Array<double, Dynamic, 1> d1 = 1 / (1 + exp(cut2)) - 1 / (1 - exp(cuts_y1 - cuts_y2));
    Array<double, Dynamic, 1> d2 = 1 / (1 - exp(cuts_y2 - cuts_y1)) - 1 / (1 + exp(cut1));

    if (!is_constant_struct<T_x>::value && !is_constant_struct<T_beta>::value) {
      Matrix<double, 1, Dynamic> deriv_common = d1 - d2;
      if (!is_constant_struct<T_x>::value) {
        ops_partials.edge1_.partials_ = (beta_val_vec * deriv_common).transpose();
      }
      if (!is_constant_struct<T_beta>::value) {
        //TODO single transpose?
        ops_partials.edge2_.partials_ = x_val.transpose() * deriv_common.transpose();
      }
    }
    if (!is_constant_struct<T_cuts>::value) {
      //TODO iter over insts?
      for (int c = 0; c < N_classes-1; c++) {
        ops_partials.edge3_.partials_[c] = y.cwiseEqual(c + 1).select(d2, 0).sum() - y.cwiseEqual(c + 2).select(d1, 0).sum();
      }
    }
  }
  return ops_partials.build(logp);
}

template<typename T_y, typename T_x, typename T_beta, typename T_cuts>
typename stan::return_type<T_x, T_beta, T_cuts>::type
ordered_logistic_glm_lpmf(
        const T_y& y, const Eigen::Matrix<T_x, Eigen::Dynamic, Eigen::Dynamic>& x,
        const T_beta& beta, const T_cuts& cuts) {
  ordered_logistic_glm_lpmf<false>(y,x,beta,cuts);
}

}
}

#endif