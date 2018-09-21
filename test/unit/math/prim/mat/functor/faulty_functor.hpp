#ifndef TEST_UNIT_MATH_PRIM_MAT_FUNCTOR_FAULTY_FUNCTOR_HPP
#define TEST_UNIT_MATH_PRIM_MAT_FUNCTOR_FAULTY_FUNCTOR_HPP

#include <vector>
#include <stdexcept>

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
    if (stan::math::abs(theta(0) + 1.0) < 1E-7) {
      throw std::domain_error("Illegal parameter!");
    }
    return (res);
  }
};

#endif
