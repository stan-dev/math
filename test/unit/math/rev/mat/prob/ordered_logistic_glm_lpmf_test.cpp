#include <stan/math.hpp>
#include <gtest/gtest.h>
#include <Eigen/Core>
#include <iostream>

//TODO cleanup
using std::cout;
using std::endl;

using stan::math::var;
using stan::math::ordered_logistic_lpmf;
using stan::math::ordered_logistic_glm_lpmf;
using Eigen::MatrixXd;
using Eigen::Matrix;
using Eigen::Array;
using Eigen::Dynamic;
using Eigen::VectorXd;
using Eigen::RowVectorXd;

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

template <bool propto, typename T_x, typename T_beta, typename T_cuts>
typename stan::return_type<T_x, T_beta, T_cuts>::type
ordered_logistic_glm_simple_lpmf(
        const Matrix<int, Dynamic, 1>& y, const Matrix<T_x, Dynamic, Dynamic>& x,
        const T_beta& beta, const T_cuts& cuts) {
  typedef typename stan::return_type<T_x, T_beta>::type T_x_beta;
  using stan::math::as_column_vector_or_scalar;

  auto& beta_col = as_column_vector_or_scalar(beta);

  Eigen::Matrix<T_x_beta, Dynamic, 1> location = x.template cast<T_x_beta>() * beta_col.template cast<T_x_beta>();

  return ordered_logistic_lpmf<propto>(y, location, cuts);
}


double ordered_logistic_lpmf2(const Matrix<int, Dynamic, 1>& y, const VectorXd& location, const VectorXd& cuts){
  using stan::math::log1p_exp;
  using stan::math::log1m_exp;
  using stan::math::log_inv_logit_diff;
  using std::abs;
  using std::log1p;

  double logp=0;
  int N = y.size();
  int C = cuts.size()+1;

  VectorXd cuts_y1(N), cuts_y2(N);
  for (int i = 0; i < N; i++) {
    if(y[i]!=C)
      cuts_y1[i] = cuts[y[i] - 1];
    if(y[i]!=1)
      cuts_y2[i] = cuts[y[i] - 2];
  }

  Array<double,Dynamic,1> cut2 = location - cuts_y2;
  Array<double,Dynamic,1> cut1 = location - cuts_y1;

  Array<double,Dynamic,1> m_log_1p_exp_cut1 = -cut1 * (cut1 > 0.0).cast<double>() - (-cut1.abs()).exp().log1p();
  Array<double,Dynamic,1> m_log_1p_exp_m_cut2 = cut2 * (cut2 <= 0.0).cast<double>() - (-cut2.abs()).exp().log1p();
  //TODO common subexpressions?
  logp = y.cwiseEqual(1).select(m_log_1p_exp_cut1,
                                y.cwiseEqual(C).select(m_log_1p_exp_m_cut2,
                                                       m_log_1p_exp_m_cut2 + log1m_exp(cut1 - cut2).array() + m_log_1p_exp_cut1
                                )).sum();

//  logp = y.cwiseEqual(1).select(-cut2 * (cut2 > 0.0).cast<double>() - (-abs(cut2)).exp().log1p(),
//         y.cwiseEqual(C).select(cut1 * (cut1 <= 0.0).cast<double>() - (-abs(cut1)).exp().log1p(),
//                                cut1 * (cut1 <= 0.0).cast<double>() - (-cut1.abs()).exp().log1p() + log1m_exp(cut2 - cut1).array() - cut2 * (cut2 > 0.0).cast<double>() - (-cut2.abs()).exp().log1p()
//                                )).sum();

//  for (int n = 0; n < N; ++n) {
////    double cut1, cut2;
////    if (y[n] == 1) {
////      cut1 = -INFINITY;
////      cut2 = cuts[0];
////    }
////    else if(y[n] == C){
////      cut1 = cuts[C - 2];
////      cut2 = INFINITY;
////    }
////    else{
////      cut1 = cuts[y[n] - 2];
////      cut2 = cuts[y[n] - 1];
////    }
////    logp += log_inv_logit_diff(location[n] - cut1,
////                               location[n] - cut2);
//    if (y[n] == 1) {
//      logp -= log1p_exp(location[n] - cuts[0]);
////      double d = inv_logit(location[n] - cuts[0]);
////
////      if (!is_constant_struct<T_loc>::value)
////        ops_partials.edge1_.partials_[n] -= d;
////
////      if (!is_constant_struct<T_cut>::value)
////        ops_partials.edge2_.partials_vec_[n](0) = d;
//
//    } else if (y[n] == C) {
//      logp -= log1p_exp(cuts[C - 2] - location[n]);
////      double d = inv_logit(cuts[C - 2] - location[n]);
////
////      if (!is_constant_struct<T_loc>::value)
////        ops_partials.edge1_.partials_[n] = d;
////
////      if (!is_constant_struct<T_cut>::value)
////        ops_partials.edge2_.partials_vec_[n](C - 2) -= d;
//
//    } else {
////      double d1
////              = inv(1 - exp(cuts[y[n] - 1] - cuts[y[n] - 2]))
////                - inv_logit(cuts[y[n] - 2] - location[n]);
////      double d2
////              = inv(1 - exp(cuts[y[n] - 2] - cuts[y[n] - 1]))
////                - inv_logit(cuts[y[n] - 1] - location[n]);
//      double a = location[n] - cuts[y[n] - 2];
//      double b = location[n] - cuts[y[n] - 1];
//      logp += a * (a <= 0.0) - log1p(exp(-abs(a))) + log1m_exp(b - a) - b * (b > 0.0) - log1p(exp(-abs(b)));
////      logp += a - log1p_exp(a) + log1m_exp(b - a) - log1p_exp(b);
////      logp += log_inv_logit_diff(a, b);
////      logp += log_inv_logit_diff(location[n] - cuts[y[n] - 2],
////                                 location[n] - cuts[y[n] - 1]);
//
////      if (!is_constant_struct<T_loc>::value)
////        ops_partials.edge1_.partials_[n] -= d1 + d2;
////
////      if (!is_constant_struct<T_cut>::value) {
////        ops_partials.edge2_.partials_vec_[n](y[n] - 2) += d1;
////        ops_partials.edge2_.partials_vec_[n](y[n] - 1) += d2;
////      }
//    }
//  }
  return logp;
}

//TEST(ProbDistributionsOrderedLogisticGLM, double_match_s){
//  double eps = 1e-13;
//  int N=5;
//  int C=3;
//  Matrix<int, Dynamic, 1> y(N);
//  y << 1,4,3,3,2;
//  VectorXd location(N);
//  location << 0.5,8,6,2,1;
//  VectorXd cuts(C);
//  cuts << 0.9,1.1,7;
//  double res1 = ordered_logistic_lpmf(y,location,cuts);
//  double res2 = ordered_logistic_lpmf2(y,location,cuts);
//  EXPECT_FLOAT_EQ(res1,res2);
//}
//TEST(ProbDistributionsOrderedLogisticGLM, double_match_big_s){
//  double eps = 1e-13;
//  int N=98;
//  int C=15;
//  Matrix<int, Dynamic, 1> y(N);
//  for(int i=0;i<N;i++){
//    y[i] = Matrix<unsigned int, Dynamic, 1>::Random(1)[0] % C + 1;
//  }
//  VectorXd location = VectorXd::Random(N);
//  VectorXd cuts = (VectorXd::Random(C).array()+1)/2/C;
//  for(int i=1;i<C;i++){
//    cuts[i]+=cuts[i-1];
//  }
//  double res1 = ordered_logistic_lpmf(y,location,cuts);
//  double res2 = ordered_logistic_lpmf2(y,location,cuts);
//  EXPECT_FLOAT_EQ(res1,res2);
//}

//TEST(ProbDistributionsOrderedLogisticGLM, glm_matches_ordered_logistic_doubles){
//  double eps = 1e-13;
//  int N=5;
//  int M=2;
//  int C=4;
//  Matrix<int, Dynamic, 1> y(N);
//  y << 1,4,3,3,2;
//  VectorXd cuts(C-1);
//  cuts << 0.9,1.1,7;
//  VectorXd beta(M);
//  beta << 1.1,0.4;
//  MatrixXd x(N,M);
//  x << 1,2,3,4,5,6,7,8,9,0;
//  EXPECT_FLOAT_EQ(ordered_logistic_glm_lpmf(y,x,beta,cuts),ordered_logistic_glm_simple_lpmf<false>(y,x,beta,cuts));
////  EXPECT_FLOAT_EQ(ordered_logistic_glm_lpmf<true>(y,x,beta,cuts),ordered_logistic_glm_simple_lpmf<true>(y,x,beta,cuts));
//}

TEST(ProbDistributionsOrderedLogisticGLM, glm_matches_ordered_logistic_vars){
  double eps = 1e-13;
  int N=5;
  int M=2;
  int C=3;
  Matrix<int, Dynamic, 1> y(N);
  y << 1,1,2,4,4;
  Matrix<var,Dynamic,1> cuts1(C), cuts2(C);
  cuts1 << 0.9,1.1,7;
  cuts2 << 0.9,1.1,7;
  Matrix<var,Dynamic,1> beta1(M), beta2(M);
  beta1 << 1.1,0.4;
  beta2 << 1.1,0.4;
  Matrix<var,Dynamic,Dynamic> x1(N,M),x2(N,M);
  x1 << 1,2,3,4,5,6,7,8,9,0;
  x2 << 1,2,3,4,5,6,7,8,9,0;
  var res1 = ordered_logistic_glm_lpmf(y,x1,beta1,cuts1);
  var res2 = ordered_logistic_glm_simple_lpmf<false>(y,x2,beta2,cuts2);
  (res1 + res2).grad();

  Matrix<double,Dynamic,Dynamic> x_adj(N,M);

  EXPECT_NEAR(res1.val(), res2.val(), eps);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      x_adj(j,i)=x2(j, i).adj();
      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
    }
  }
  for (int i = 0; i < M; i++) {
    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
  }
  for (int i = 0; i < C; i++) {
    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
  }
//  var res2 = ordered_logistic_glm_lpmf<true>(y,x,beta,cuts),ordered_logistic_glm_simple_lpmf<true>(y,x,beta,cuts);

}

//TEST(ProbDistributionsOrderedLogisticGLM, glm_matches_ordered_logistic_vars_big){
//  double eps = 1e-7;
//  int N=890;
//  int M=15;
//  int C=10;
//  Matrix<int, Dynamic, 1> y(N);
//  for (int i = 0; i < N; i++) {
//    y[i] = Matrix<unsigned int, Dynamic, 1>::Random(1)[0] % (C+1) + 1;
//  }
//  VectorXd cuts_double = (VectorXd::Random(C).array()+1)/2/C;
//  for(int i=1;i<C;i++){
//    cuts_double[i]+=cuts_double[i-1];
//  }
//  Matrix<var,Dynamic,1> cuts1 = cuts_double, cuts2 = cuts_double;
//  VectorXd beta_double = VectorXd::Random(M);
//  Matrix<var,Dynamic,1> beta1 = beta_double, beta2 = beta_double;
//  MatrixXd x_double = MatrixXd::Random(N, M);
//  Matrix<var,Dynamic,Dynamic> x1 = x_double, x2 = x_double;
//  var res1 = ordered_logistic_glm_lpmf(y,x1,beta1,cuts1);
//  var res2 = ordered_logistic_glm_simple_lpmf<false>(y,x2,beta2,cuts2);
//  (res1 + res2).grad();
//
//  Matrix<double,Dynamic,Dynamic> x_adj(N,M);
//
//  EXPECT_NEAR(res1.val(), res2.val(), eps);
//  for (int i = 0; i < M; i++) {
//    for (int j = 0; j < N; j++) {
//      x_adj(j,i)=x2(j, i).adj();
//      EXPECT_NEAR(x1(j, i).adj(), x2(j, i).adj(), eps);
//    }
//  }
//  for (int i = 0; i < M; i++) {
//    EXPECT_NEAR(beta1[i].adj(), beta2[i].adj(), eps);
//  }
//  for (int i = 0; i < C; i++) {
//    EXPECT_NEAR(cuts1[i].adj(), cuts2[i].adj(), eps);
//  }
////  var res2 = ordered_logistic_glm_lpmf<true>(y,x,beta,cuts),ordered_logistic_glm_simple_lpmf<true>(y,x,beta,cuts);
//}