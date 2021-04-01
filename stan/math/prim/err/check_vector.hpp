#ifndef STAN_MATH_PRIM_ERR_CHECK_VECTOR_HPP
#define STAN_MATH_PRIM_ERR_CHECK_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <sstream>
#include <string>
#include <typeinfo>

#ifdef STAN_OPENCL
#include <stan/math/opencl/value_type.hpp>
#endif

namespace stan {
namespace math {

/**
 * Check the input is either a row vector or column vector or
 *  a matrix with a single row or column.
 *
 * @tparam Mat Input type
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param x Input
 * @throw <code>std::invalid_argument</code> if x is not a row or column
 *   vector.
 */
template <typename Mat,
          require_any_t<is_matrix<Mat>,
                        is_prim_or_rev_kernel_expression<Mat>>* = nullptr>
inline void check_vector(const char* function, const char* name, const Mat& x) {
  if (!(x.rows() == 1 || x.cols() == 1)) {
    STAN_NO_RANGE_CHECKS_RETURN;
    [&]() STAN_COLD_PATH {
      std::ostringstream msg;
      msg << ") has " << x.rows() << " rows and " << x.cols()
          << " columns but it should be a vector so it should "
          << "either have 1 row or 1 column";
      std::string msg_str(msg.str());
      invalid_argument(function, name, 0.0, "(", msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif
