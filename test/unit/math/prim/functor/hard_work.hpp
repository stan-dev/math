#ifndef TEST_UNIT_MATH_PRIM_FUNCTOR_HARD_WORK_HPP
#define TEST_UNIT_MATH_PRIM_FUNCTOR_HARD_WORK_HPP

#include <vector>

struct hard_work {
  hard_work() {}
  template <typename T1, typename T2>
  Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> operator()(
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& eta,
      const Eigen::Matrix<T2, Eigen::Dynamic, 1>& theta,
      const std::vector<double>& x_r, const std::vector<int>& x_i,
      std::ostream* msgs = 0) const {
    using result_type = stan::return_type_t<T1, T2>;
    Eigen::Matrix<result_type, Eigen::Dynamic, 1> res;
    res.resize(2 + x_i[0]);
    res.setZero();
    res(0) = theta(0) * theta(0) + eta(0);
    res(1) = x_r[0] * theta(1) * theta(0) + 2 * eta(0) + eta(1);
    return (res);
  }
};

#endif
