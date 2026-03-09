#ifndef STAN_MATH_PRIM_FUN_BUILD_MATRIX
#define STAN_MATH_PRIM_FUN_BUILD_MATRIX

namespace stan::math {


template <typename F, typename... Args, require_all_not_st_var<Args...>* = nullptr>
auto build_matrix(int rows, int cols, F&& f, Args&&... args) {
  Eigen::MatrixXd m(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      m(i, j) = f(i, j, args...);
    }
  }
  return m;
}
}

#endif