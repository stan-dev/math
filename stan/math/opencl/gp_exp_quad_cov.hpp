#ifndef STAN_MATH_OPENCL_GP_EXP_QUAD_COV_HPP
#define STAN_MATH_OPENCL_GP_EXP_QUAD_COV_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernels/gp_exp_quad_cov.hpp>
#include <stan/math/opencl/err/check_matching_dims.hpp>
#include <CL/cl.hpp>

namespace stan {
namespace math {
/**
 * Squared exponential kernel on the GPU.
 *
 * @param x input vector
 * @param sigma sigma
 * @param length_scale
 *
 * @return Squared exponential kernel of x.
 */
inline matrix_cl gp_exp_quad_cov(matrix_cl& x, double sigma, double length_scale) {
  matrix_cl res(x.cols(), x.cols());
  try {
    opencl_kernels::gp_exp_quad_cov(cl::NDRange(x.cols(), x.cols()),
                                    x.buffer(), res.buffer(),
                                    sigma * sigma, -0.5 / (length_scale * length_scale), x.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exp_quad_cov", e);
  }
  return res;
}

/**
 * Squared exponential kernel on the GPU.
 *
 * @param x first input vector
 * @param y second input vector
 * @param sigma sigma
 * @param length_scale
 *
 * @return Squared exponential kernel of x.
 */
inline matrix_cl gp_exp_quad_cov(matrix_cl& x, matrix_cl& y, double sigma, double length_scale) {
  check_size_match("gp_exp_quad_cov", "x", x.rows(), "y", y.rows());
  matrix_cl res(x.cols(), y.cols());
  try {
    opencl_kernels::gp_exp_quad_cov_2(cl::NDRange(x.cols(), y.cols()),
                                      x.buffer(), y.buffer(), res.buffer(),
                                      sigma * sigma, -0.5 / (length_scale * length_scale), x.cols(), y.cols(), x.rows());
  } catch (const cl::Error& e) {
    check_opencl_error("gp_exp_quad_cov_2", e);
  }
  return res;
}

/*
inline
typename Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
        Eigen::Dynamic, Eigen::Dynamic>
gp_exp_quad_cov(const std::vector<T_x>& x, const T_sigma& sigma,
                const T_l& length_scale) {
    using std::exp;
    check_positive("gp_exp_quad_cov", "marginal variance", sigma);
    check_positive("gp_exp_quad_cov", "length-scale", length_scale);
    for (size_t n = 0; n < x.size(); ++n)
        check_not_nan("gp_exp_quad_cov", "x", x[n]);

    Eigen::Matrix<typename stan::return_type<T_x, T_sigma, T_l>::type,
            Eigen::Dynamic, Eigen::Dynamic>
            cov(x.size(), x.size());

    int x_size = x.size();
    if (x_size == 0)
        return cov;

    const auto& x_val = value_of(x);
    const auto sigma_val = value_of(sigma);
    const auto length_scale_val = value_of(length_scale);

    auto sigma_sq = square(sigma_val);
    auto neg_half_inv_l_sq = -0.5 / square(length_scale_val);

    operands_and_partials<std::vector<T_x>, T_sigma, T_l> ops_partials(x, sigma, length_scale);

    for (int j = 0; j < (x_size - 1); ++j) {
        cov(j, j) = sigma_sq;
        for (int i = j + 1; i < x_size; ++i) {
            double distance = x_val[i] - x_val[j];
            double sq_distance = distance * distance;
            cov(i, j) = sigma_sq * exp(sq_distance * neg_half_inv_l_sq);
            cov(j, i) = cov(i, j);
            if (!is_constant_struct<T_x>::value) {
                ops_partials.edge1_.partials_ = ;
            }
        }
    }
    cov(x_size - 1, x_size - 1) = sigma_sq;
    return cov;
}*/

}  // namespace math
}  // namespace stan

#endif
#endif
