#pragma once

#include <stdexcept>
#include <vector>

struct faulty_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<typename stan::return_type<T1, T2>::type, Eigen::Dynamic, 1>
  operator()(const Eigen::Matrix<T1, Eigen::Dynamic, 1>& eta,
             const Eigen::Matrix<T2, Eigen::Dynamic, 1>& theta,
             const std::vector<double>& x_r, const std::vector<int>& x_i,
             std::ostream* msgs = 0) const {
    typedef typename stan::return_type<T1, T2>::type result_type;
    Eigen::Matrix<result_type, Eigen::Dynamic, 1> res;
    res.resize(2);
    res(0) = theta(0) * theta(0);
    res(1) = x_r[0] * theta(1) * theta(0);
    // std::cout << "theta(0) = " << theta(0) << std::endl;
    if (stan::math::abs(theta(0) + 1.0) < 1E-7) {  // check for param being 1.0
      // std::cout << "THROWING! DEAL WITH IT!" << std::endl;
      throw std::domain_error("Illegal parameter!");
    }
    return (res);
  }
};
