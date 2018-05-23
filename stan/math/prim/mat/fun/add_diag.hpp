#ifndef STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP
#define STAN_MATH_PRIM_MAT_FUN_ADD_DIAG_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/fun/fmin.hpp>

// testing
#include </home/andre/math/stan/math/prim/scal/meta/return_type.hpp>

namespace stan {
namespace math {

template <typename T_m, typename T_a>
inline typename Eigen::Matrix<typename return_type<T_m, T_a>::type, -1, -1>
add_diag(const Eigen::Matrix<T_m, -1, -1> &mat, const T_a &to_add) {
  // check matrix non-empty
  // check nrow, ncol >= 1
  // check to_add not nan, not inf, 

  size_t length_diag;
  length_diag = fmin(mat.rows(), mat.cols()); // or use Eigen's function?
  //  Eigen::Matrix<typename return_type<T_m, T_a>::type, -1, -1> out;
  Eigen::MatrixXd out;
  //  out = mat; //  I could likewise remote const and motify mat itself
  //  Eigen::Matrix<typename return_type<T_m, T_a>::type, 1, -1> temp_vec;

  // Eigen::Matrix<T_m, 1, -1> temp_vec;
  // temp_vec = mat.diagonal(); //  so we need not call .diagonal() more than once
  
  // for (size_t i = 0; i < length_diag; ++i)
  //   out(i, i) = temp_vec[i] + to_add;
    //    temp_vec[i] += to_add;
  //  return out;
  return mat;
}

}
}


#endif
