#pragma once

struct hard_work {
  template<typename T1, typename T2>
  Eigen::Matrix<typename stan::return_type<T1, T2>::type, Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& eta, const Eigen::Matrix<T2, Eigen::Dynamic, 1>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i) const {
    typedef typename stan::return_type<T1, T2>::type result_type;
    Eigen::Matrix<result_type, Eigen::Dynamic, 1> res(2);
    res(0) = theta(0)*theta(0);
    res(1) = x_r[0]*theta(1)*theta(0);
    return(res);
  }

  template<typename T1, typename T2>
  static
  Eigen::Matrix<typename stan::return_type<T1, T2>::type, Eigen::Dynamic, 1>
  apply(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& eta, const Eigen::Matrix<T2, Eigen::Dynamic, 1>& theta, const std::vector<double>& x_r, const std::vector<int>& x_i) {
    const hard_work f;
    return f(eta, theta, x_r, x_i);
  }
};
