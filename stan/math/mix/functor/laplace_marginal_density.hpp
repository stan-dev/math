#ifndef STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#define STAN_MATH_MIX_FUNCTOR_LAPLACE_MARGINAL_DENSITY_HPP
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/mix/functor/laplace_likelihood.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/rev/functor.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/quad_form_diag.hpp>
#include <stan/math/prim/functor/iter_tuple_nested.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>
#include <optional>
//#define LAPLACE_DEBUG
#ifdef LAPLACE_DEBUG
#include <iomanip>
#endif
/**
 * @file
 * Reference for calculations of marginal and its gradients:
 * Margossian et al (2020), https://arxiv.org/abs/2004.12550
 * and Margossian (2023), https://arxiv.org/pdf/2306.14976
 */

namespace stan {
namespace math {

struct laplace_line_search_options {
  double c1{1e-4};
  double c2{0.9};
  double tau{0.5};
  double min_alpha{1e-12};
  double abs_grad_threshold{1e-14};
  double abs_obj_threshold{1e-14};
};

/**
 * Options for the laplace sampler
 */
struct laplace_options_base {
  /* Size of the blocks in block diagonal hessian*/
  int hessian_block_size{1};
  /**
   * Which Newton solver to use: (B matrix in equation 1 of
   * https://arxiv.org/pdf/2306.14976) (1) method using the cholesky
   * decomposition of `W` (the negative Hessian of log likelihood) (2) method
   * using the cholesky decomposition of `K` (the covariance matrix) (3) method
   * using an LU decomposition (more general, but slower)
   */
  int solver{1};
  /* Maximum number of steps in line search */
  int max_steps_line_search{100};
  /* iterations end when difference in objective function is less than tolerance
   */
  double tolerance{1e-12};
  /* Maximum number of steps*/
  int max_num_steps{100};
  laplace_line_search_options line_search;
};

template <bool HasInitTheta>
struct laplace_options;

template <>
struct laplace_options<false> : public laplace_options_base {
};

template <>
struct laplace_options<true> : public laplace_options_base {
  /* Value for user supplied initial theta  */
  Eigen::VectorXd theta_0{0};
};

using laplace_options_default = laplace_options<false>;
using laplace_options_user_supplied = laplace_options<true>;
namespace internal {

template <typename Covar, typename ThetaVec, typename WR, typename L_t,
          typename A_vec, typename ThetaGrad, typename LU_t, typename KRoot>
struct laplace_density_estimates {
  /* log marginal density */
  double lmd{std::numeric_limits<double>::infinity()};
  /* Evaluated covariance function for the latent gaussian variable */
  Covar covariance;
  /* ThetaVec at the mode */
  ThetaVec theta;
  /* negative hessian or sqrt of negative hessian */
  WR W_r;
  /* Lower left of cholesky decomposition of stabilized inverse covariance */
  L_t L;
  /* inverse covariance times theta at the mode */
  A_vec a;
  /* the gradient of the log density with respect to theta */
  ThetaGrad theta_grad;
  /* LU matrix from solver 3 */
  LU_t LU;
  /* Cholesky of the covariance matrix */
  KRoot K_root;
  laplace_density_estimates(double lmd_, Covar&& covariance_, ThetaVec&& theta_,
                            WR&& W_r_, L_t&& L_, A_vec&& a_,
                            ThetaGrad&& theta_grad_, LU_t&& LU_,
                            KRoot&& K_root_)
      : lmd(lmd_),
        covariance(std::move(covariance_)),
        theta(std::move(theta_)),
        W_r(std::move(W_r_)),
        L(std::move(L_)),
        a(std::move(a_)),
        theta_grad(std::move(theta_grad_)),
        LU(std::move(LU_)),
        K_root(std::move(K_root_)) {}
};

/**
 * Returns the principal square root of a block diagonal matrix.
 * @tparam WRootMat A type inheriting from `Eigen::EigenBase`.
 * @param W_root The output matrix to store the square root.
 * @param W The input block diagonal matrix.
 * @param block_size The size of each block in the block diagonal matrix.
 */
template <typename WRootMat>
inline void block_matrix_sqrt(WRootMat& W_root,
                              const Eigen::SparseMatrix<double>& W,
                              const Eigen::Index block_size) {
  int n_block = W.cols() / block_size;
  Eigen::MatrixXd local_block(block_size, block_size);
  Eigen::MatrixXd local_block_sqrt(block_size, block_size);
  Eigen::MatrixXd sqrt_t_mat = Eigen::MatrixXd::Zero(block_size, block_size);
  // No block operation available for sparse matrices, so we have to loop
  // See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title7
  for (int i = 0; i < n_block; i++) {
    sqrt_t_mat.setZero();
    local_block
        = W.block(i * block_size, i * block_size, block_size, block_size);
    if (Eigen::isnan(local_block.array()).any()) {
      throw std::domain_error(
          std::string("Error in block_matrix_sqrt: "
                      "NaNs detected in block diagonal starting at (")
          + std::to_string(i) + ", " + std::to_string(i) + ")");
    }
    // Issue here, sqrt is done over T of the complex schur
    Eigen::RealSchur<Eigen::MatrixXd> schurOfA(local_block);
    // Compute Schur decomposition of arg
    const auto& t_mat = schurOfA.matrixT();
    const auto& u_mat = schurOfA.matrixU();
    // Check if diagonal of schur is not positive
    if ((t_mat.diagonal().array() < 0).any()) {
      throw std::domain_error(
          std::string("Error in block_matrix_sqrt: "
                      "values less than 0 detected in block diagonal's schur "
                      "decomposition starting at (")
          + std::to_string(i) + ", " + std::to_string(i) + ")");
    }
    try {
      // Compute square root of T
      Eigen::matrix_sqrt_quasi_triangular(t_mat, sqrt_t_mat);
      // Compute square root of arg
      local_block_sqrt = u_mat * sqrt_t_mat * u_mat.adjoint();
    } catch (const std::exception& e) {
      throw std::domain_error(
          "Error in block_matrix_sqrt: "
          "The matrix is not positive definite");
    }
    for (int k = 0; k < block_size; k++) {
      for (int j = 0; j < block_size; j++) {
        W_root.coeffRef(i * block_size + j, i * block_size + k)
            = local_block_sqrt(j, k);
      }
    }
  }
}

/**
 * @brief Performs a Cholesky decomposition on a block diagonal matrix.
 * @tparam WRootMat A type inheriting from `Eigen::EigenBase`.
 * @param W_root The output matrix to store the square root.
 * @param W The input block diagonal matrix.
 * @param block_size The size of each block in the block diagonal matrix.
 */
template <typename WRootMat>
inline void block_matrix_chol_L(WRootMat& W_root,
                                const Eigen::SparseMatrix<double>& W,
                                const Eigen::Index block_size) {
  int n_block = W.cols() / block_size;
  Eigen::MatrixXd local_block(block_size, block_size);
  Eigen::MatrixXd local_block_sqrt(block_size, block_size);
  Eigen::MatrixXd sqrt_t_mat = Eigen::MatrixXd::Zero(block_size, block_size);
  // No block operation available for sparse matrices, so we have to loop
  // See https://eigen.tuxfamily.org/dox/group__TutorialSparse.html#title7
  for (int i = 0; i < n_block; i++) {
    sqrt_t_mat.setZero();
    local_block
        = W.block(i * block_size, i * block_size, block_size, block_size);
    if (Eigen::isnan(local_block.array()).any()) {
      throw std::domain_error(
          std::string("Error in block_matrix_sqrt: "
                      "NaNs detected in block diagonal starting at (")
          + std::to_string(i) + ", " + std::to_string(i) + ")");
    }
    try {
      // Compute square root of T
      Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt(local_block);
      if (llt.info() != Eigen::Success) {
        throw std::runtime_error("Cholesky failed on block "
                                 + std::to_string(i));
      }
      const auto Lb = llt.matrixL();
      for (int k = 0; k < block_size; k++) {
        for (int j = k; j < block_size; j++) {
          W_root.coeffRef(i * block_size + j, i * block_size + k) = Lb(j, k);
        }
      }
    } catch (const std::exception& e) {
      // As a backup do the schur decomposition for this block diagonal
      local_block
          = W.block(i * block_size, i * block_size, block_size, block_size);
      // Issue here, sqrt is done over T of the complex schur
      Eigen::RealSchur<Eigen::MatrixXd> schurOfA(local_block);
      // Compute Schur decomposition of arg
      const auto& t_mat = schurOfA.matrixT();
      const auto& u_mat = schurOfA.matrixU();
      // Check if diagonal of schur is not positive
      if ((t_mat.diagonal().array() < 0).any()) {
        throw std::domain_error(
            std::string("Error in block_matrix_sqrt: "
                        "values less than 0 detected in block diagonal's schur "
                        "decomposition starting at (")
            + std::to_string(i) + ", " + std::to_string(i) + ")");
      }
      try {
        // Compute square root of T
        Eigen::matrix_sqrt_quasi_triangular(t_mat, sqrt_t_mat);
        // Compute square root of arg
        local_block_sqrt.noalias() = u_mat * sqrt_t_mat * u_mat.adjoint();
      } catch (const std::exception& e) {
        throw std::domain_error(
            "Error in block_matrix_sqrt: "
            "The matrix is not positive definite");
      }
      for (int k = 0; k < block_size; k++) {
        for (int j = 0; j < block_size; j++) {
          W_root.coeffRef(i * block_size + j, i * block_size + k)
              = local_block_sqrt(j, k);
        }
      }
      throw std::domain_error(
          "Error in block_matrix_sqrt: "
          "The matrix is not positive definite");
    }
  }
}

/**
 * Collect the adjoints from the input and add them to the output.
 * @tparam ZeroInput If true, the adjoints of the input will be set to zero
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `var` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <bool ZeroInput = false, typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_all_var_scalar<Input>>* = nullptr>
inline void collect_adjoints(Output& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<var, -1, 1>> input_map(input_i.data(),
                                                          input_i.size());
          output_map.array() += input_map.adj().array();
          if constexpr (ZeroInput) {
            input_map.adj().setZero();
          }
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += input_i.adj().array();
          if constexpr (ZeroInput) {
            input_i.adj().setZero();
          }
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += input_i.adj();
          if constexpr (ZeroInput) {
            input_i.adj() = 0;
          }
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}

/**
 * Throws an error if the parameter contains NaN or Inf values.
 * @tparam NameStr Type of the name string, e.g. `std::string` or `char*`.
 * @tparam ParamStr Type of the parameter string, e.g. `std::string` or `char*`.
 * @tparam Param Type of the parameter such as a vector, matrix, or scalar.
 * @param name_str Name of the function or context where the error occurred.
 * @param param_str Name of the parameter that contains NaN or Inf values.
 * @param param The parameter to check for NaN or Inf values.
 */
template <typename NameStr, typename ParamStr, typename Param>
inline STAN_COLD_PATH void throw_nan(NameStr&& name_str, ParamStr&& param_str,
                                     Param&& param) {
  std::string msg = std::string("Error in ") + name_str + ": "
                    + std::string(param_str) + " contains NaN values";
  if ((Eigen::isnan(param.array()) || Eigen::isinf(param.array())).all()) {
    msg += " for all values.";
    throw std::domain_error(msg);
  }
  msg += " at indices [";
  for (int i = 0; i < param.size(); ++i) {
    if (std::isnan(param(i) || std::isinf(param(i)))) {
      msg += std::to_string(i) + ", ";
    }
  }
  msg.pop_back();
  msg.pop_back();
  msg += "].";
  throw std::domain_error(msg);
}

/**
 * @brief  One–dimensional cubic interpolation helper used by the Wolfe
 *         zoom phase to pick the step length that *maximises* an
 *         auxiliary cubic model *g* on the interval \f$\,0 < x <
 * \text{hi\_limit}\f$.
 *
 * Four pieces of data are provided:
 *
 *   * the slope at the left end (`dir_lo  = g'(0)`);
 *   * the slope at the right end (`dir_high = g'(hi_bound)`);
 *   * the function value difference `obj_diff = g(hi_bound) − g(0)`; and
 *   * the abscissa of the right end (`hi_bound`).
 *
 * A cubic polynomial that matches those values is constructed analytically
 *     \f[
 *       g'(x) = a_3 x^2 + a_2 x + a_1 ,
 *     \f]
 * and its two stationary points *inside* the current bracket plus the
 * upper bound itself are tested.  The point that gives the largest
 * \f$g(x)\f$ is returned.
 *
 * **Numerical safeguards:**
 * * **Degenerate cubic → midpoint.**  If the leading coefficient
 *   \f$a_3\f$ is almost zero (`|a3| < 1 e−14`) the polynomial becomes
 *   quadratic or linear; the function falls back to the simpler and
 *   numerically safer choice `hi_limit / 2`.
 * * **Negative discriminant clamp.**  Round‑off may make the quantity
 *   under the square‑root slightly negative; it is clamped to zero so the
 *   code never feeds NaNs into `std::sqrt`.
 * * **Safety clipping.**  Only roots that lie strictly inside
 *   \f$(0,\text{hi\_limit})\f$ compete with `hi_limit` itself.  This guards
 *   against wandering outside the current Wolfe bracket.
 *
 * @tparam T A scalar type
 * @param[in] dir_lo   \f$g'(0)\f$ – directional derivative at the **left** end
 * of the current interval.
 * @param[in] hi_bound Length of the interval, \f$\text{hi\_bound}=
 * \alpha_{\text{high}} - \alpha_{\text{low}}\f$.
 * @param[in] obj_diff Difference in the *primitive* between the two anchors,
 * \f$g(\text{hi\_bound}) - g(0)\f$.
 * @param[in] dir_high \f$g'(\text{hi\_bound})\f$ – directional derivative at
 * the **right** end of the interval.
 * @param[in] hi_limit upper limit on the step that may be returned
 *
 * @return The abscissa \f$x^{*}\in(0,\text{hi\_limit})\f$ where the cubic
 *         predicts the largest ascent of *g*.  If numerical safeguards
 *         are triggered, the midpoint `hi_limit / 2` is returned.
 *       loops.
 */
template <typename T>
inline auto cubic_interp_max(T dir_lo, T hi_bound, T obj_diff, T dir_high,
                             T hi_limit) {
  // g′(x) = a₃x² + a₂x + a₁  (with g(0)=0)
  const T a3 = (-12.0 * obj_diff + 6.0 * hi_bound * (dir_lo + dir_high))
               / std::pow(hi_bound, 3);
  const T a2 = 6.0 * obj_diff / (hi_bound * hi_bound)
               - (4.0 * dir_lo + 2.0 * dir_high) / hi_bound;
  const T a1 = dir_lo;
  // Degenerate or nearly-linear case → midpoint
  if (std::abs(a3) < 1e-14) {
    return hi_limit * T(0.5);
  }
  T disc2 = a2 * a2 - 2.0 * a1 * a3;
  disc2 = (disc2 < 0) ? 0 : disc2;
  const T disc = std::sqrt(disc2);
  const T r1 = (-a2 + disc) / a3;
  const T r2 = (-a2 - disc) / a3;
  auto g = [&](T x) { return ((a3 / 3) * x + (a2 / 2)) * x * x + a1 * x; };
  T best_x = 0, best_val = g(0);
  for (T r : {r1, r2, hi_limit}) {
    if (r > 0 && r < hi_limit) {
      T v = g(r);
      const bool v_best = v > best_val;
      best_x = v_best ? r : best_x;
      best_val = v_best ? v : best_val;
    }
  }
  return best_x;
}

namespace debug {

constexpr void print(int tabs) {}

#ifdef LAPLACE_DEBUG
template <typename Val, typename... Types>
void print(int tabs, const char* name, Val&& val, Types&&... args) {
  for (int i = 0; i < tabs; ++i) {
    std::cout << "\t";
  }
  std::cout << name << ": " << std::setprecision(13) << std::fixed << stan::math::eval(val) << "\n";
  print(tabs, std::forward<Types>(args)...);
}

template <typename Val, typename... Types>
void print(const char* msg, int tabs, const char* name, Val&& val, Types&&... args) {
  std::cout << msg << "\n";
  for (int i = 0; i < tabs; ++i) {
    std::cout << "\t";
  }
  std::cout << name << ": " << std::setprecision(13) << std::fixed << stan::math::eval(val) << "\n";
  print(tabs, std::forward<Types>(args)...);
}

void print(const char* msg) {
  std::cout << msg << "\n";
}

template <typename Val>
void print(const char* msg, Val&& val) {
  std::cout << msg << std::setprecision(13) << std::fixed << stan::math::eval(val) << "\n";
}


#else
template <typename Val, typename... Types>
constexpr void print(int tabs, const char* name, Val&& val, Types&&... args) {}

template <typename Val, typename... Types>
constexpr void print(const char* msg, int tabs, const char* name, Val&& val, Types&&... args) {}

constexpr void print(const char* msg) {}

template <typename Val>
constexpr void print(const char* msg, Val&& val) {}
#endif
}

template <typename Option>
inline auto check_armijo(double obj_next, double obj_init, double alpha_next,
                         double dir0, Option&& opt) {
  debug::print("\n\tcheck_armijo: ", 2,
               "success:    ", (obj_next >= obj_init + alpha_next * dir0 * opt.line_search.c1 ? "true" : "false"),
               "obj_next:   ", obj_next,
               "obj_init:   ", obj_init,
               "alpha_next: ", alpha_next,
               "dir0:       ", dir0,
               "c1:         ", opt.line_search.c1,
               "obj + alpha * dir0 * c1: ", (obj_init + alpha_next * dir0 * opt.line_search.c1)
               );
  return obj_next >= obj_init + alpha_next * dir0 * opt.line_search.c1;
}

template <typename Option>
inline auto check_wolfe_curve(double dir_deriv_next, double dir_deriv_init,
                              Option&& opt) {
  debug::print("check_wolfe_curve: ", 2,
               "success:    ", (std::abs(dir_deriv_next) <= (opt.line_search.c2 * std::abs(dir_deriv_init)) ? "true" : "false"),
               "deriv_next: ", dir_deriv_next,
               "deriv_init: ", dir_deriv_init,
               "c2:         ", opt.line_search.c2,
               "abs(d_next):   ", std::abs(dir_deriv_next),
               "abs(d_init)*c2 ", (std::abs(dir_deriv_init) * opt.line_search.c2)
              );
  return std::abs(dir_deriv_next) <= (opt.line_search.c2 * std::abs(dir_deriv_init));
}

template <typename Scalar>
Scalar CubicInterp(const Scalar &df0, const Scalar &x1, const Scalar &f1,
                   const Scalar &df1, const Scalar &loX, const Scalar &hiX) {
  const Scalar c3((-12 * f1 + 6 * x1 * (df0 + df1)) / (x1 * x1 * x1));
  const Scalar c2(-(4 * df0 + 2 * df1) / x1 + 6 * f1 / (x1 * x1));
  const Scalar &c1(df0);

  const Scalar t_s = std::sqrt(c2 * c2 - 2.0 * c1 * c3);
  const Scalar s1 = -(c2 + t_s) / c3;
  const Scalar s2 = -(c2 - t_s) / c3;

  Scalar tmpF;
  Scalar minF, minX;

  // Check value at lower bound
  minF = loX * (loX * (loX * c3 / 3.0 + c2) / 2.0 + c1);
  minX = loX;

  // Check value at upper bound
  tmpF = hiX * (hiX * (hiX * c3 / 3.0 + c2) / 2.0 + c1);
  if (tmpF < minF) {
    minF = tmpF;
    minX = hiX;
  }

  // Check value of first root
  if (loX < s1 && s1 < hiX) {
    tmpF = s1 * (s1 * (s1 * c3 / 3.0 + c2) / 2.0 + c1);
    if (tmpF < minF) {
      minF = tmpF;
      minX = s1;
    }
  }

  // Check value of second root
  if (loX < s2 && s2 < hiX) {
    tmpF = s2 * (s2 * (s2 * c3 / 3.0 + c2) / 2.0 + c1);
    if (tmpF < minF) {
      minF = tmpF;
      minX = s2;
    }
  }

  return minX;
}

template <typename Scalar>
Scalar CubicInterp(const Scalar &x0, const Scalar &f0, const Scalar &df0,
                   const Scalar &x1, const Scalar &f1, const Scalar &df1,
                   const Scalar &loX, const Scalar &hiX) {
  return x0 + CubicInterp(df0, x1 - x0, f1 - f0, df1, loX - x0, hiX - x0);
}

enum class wolfe_return {
  PASS,
  ARMIJO_PASS,
  FAIL,
  GRAD_CONV
};

/**
 * @brief  Strong‑Wolfe line search with cubic‑interpolation "zoom" for Laplace‐style log‑likelihood problems.
 *
 * This routine searches along the space of the latent gaussian *a*
 * \f$a(\alpha) = a_prev + \alpha p`, `p = a − a_prev\f$,
 * looking for the largest step \f$\alpha\f$ that satisfies the **strong‑Wolfe**
 * conditions
 *
 * \f{align*}{
 *   \phi(\alpha) &\;\ge\; \phi(0) + c_1 \alpha \phi'(0) \quad\text{(Armijo)},\\
 *   |\phi'(\alpha)| &\;\le\; c_2 |\phi'(0)| \quad\text{(curvature)}, \f}
 *
 * where
 *  \f$\phi(\alpha)=\text{obj\_fun}\bigl(a(\alpha),\;\theta(\alpha)\bigr)\f$,
 *  \f$\theta(\alpha)=\text{covariance}\;·\;a(\alpha)\f$,
 *  \f$\phi'(\alpha)=\nabla\phi(\alpha)^{\!T} p\f$.
 *
 * The search proceeds in three phases
 *
 *  1. **Back‑tracking** – halve the initial \f$\alpha\f$ until both Wolfe conditions pass.
 *  2. If the \f$alpha\f$ is 1 or goes below `min_alpha`, then we end the search early, as Laplace problems commonly accept a full Newton step.
 *  3. **Bracketing by doubling** – starting from that good \f$\alpha\f$, double the step until Armijo fails; the last good point is the left end of the bracket, the first failing point the right end.
 *  4. **Cubic zoom** – repeatedly fit a cubic through the bracket end‑points
 *     (`cubic_interp_max`), evaluate the objective/gradient at the
 *     predicted maximiser, and shrink the bracket until a Wolfe‑compliant
 *     step is found or the interval width falls below `opt.line_search.min_alpha`.
 *
 * * **Gradient reuse** – the caller provides `theta_grad` computed at the
 *   starting point; inside the loop it is overwritten with fresh gradients
 *   via `laplace_likelihood::theta_grad`.
 * * **Early‑exit for \f$\alpha\f$ = 1 or < `min_alpha`** Laplace problems commonly accept a full Newton step. The function short‑circuits to avoid any extra work.
 *
 * @tparam F Callable type of the raw log‑likelihood (passed to `laplace_likelihood::theta_grad`)
 * @tparam Obj Callable returning the scalar objective value `obj_fun(a, theta)`
 * @tparam Grad Callable returning the gradient in *a‑space* given `(a, theta, theta_grad)`
 * @tparam LLArgs Struct or tuple holding additional arguments forwarded to `ll_fun` and `theta_grad`
 * @tparam Stream Any type that implements the stream interface used by `laplace_likelihood` for diagnostic messages; may be `std::ostream`‐like or `nullptr`
 * @tparam Options Struct holding search parameters (`c1`, `c2`, `min_alpha`, ...)
 *
 * @param[in,out] theta Current \f$\theta=\Sigma a\f$; overwritten with the accepted value on success
 * @param[in,out] obj_init Objective at \f$\alpha\f$ = 0 on entry; updated to the objective at the accepted step on exit
 * @param[in,out] alpha_init Step size suggestion on input; on success holds the step that satisfied strong‑Wolfe
 * @param[in,out] theta_grad Gradient of the log‑likelihood wrt theta at the **current** theta.  Recomputed internally and returned at the final theta
 * @param[in,out] a Working copy of *a*; on success contains the accepted point
 * @param[in] a_prev Starting point (\f$\alpha\f$ = 0)
 * @param[in] ll_fun Raw log‑likelihood functor
 * @param[in] obj_fun Objective functor to be maximised
 * @param[in] grad_fun Functor returning gradint of objective with respect to `a`.
 * @param[in] covariance Symmetric positive‑definite $\Sigma$ converting *a* to theta
 * @param[in] ll_args Extra arguments forwarded to `ll_fun` and `theta_grad`
 * @param[in] opt Line‑search constants (`c1`, `c2`, etc.)
 * @param[in,out] msgs Optional diagnostics stream.  May be `nullptr` to suppress messages
 *
 * @return `true`  if a step satisfying both Wolfe conditions was found
 *         `false` if only Armijo is satisfied (strong‑Wolfe failed but a
 *                “least‑bad’’ step was returned).
 *
 * @warning The helper `cubic_interp_max` assumes its first point
 *          corresponds to *g*(0)=0; do **not** alter the baseline
 *          initialisation logic or the interpolation will become
 *          inconsistent.
 */
template <typename F, class Obj, class Grad, typename LLArgs, typename Stream,
          typename Options>
inline wolfe_return wolfe_line_search(Eigen::VectorXd& theta, double& obj_init,
                              double& alpha_init, Eigen::VectorXd& theta_grad,
                              Eigen::VectorXd& a, const Eigen::VectorXd& a_prev,
                              F&& ll_fun, Obj&& obj_fun, Grad&& grad_fun,
                              const Eigen::MatrixXd& covariance,
                              LLArgs&& ll_args, Options&& opt, Stream* msgs) {
  Eigen::VectorXd p = a - a_prev;
  const Eigen::VectorXd g0 = -covariance * a_prev + covariance * theta_grad;
  double dir_deriv_init = g0.dot(p);
  auto update_step_vals = [&p, &a_prev, &covariance,
    &ll_fun, &ll_args, &obj_fun, &grad_fun, msgs](
    auto& a_in, auto& theta_in, auto& theta_grad_in,
    auto& obj_in, auto& dir_deriv_in, auto& alpha) {
      a_in        = a_prev + alpha * p;
      theta_in    = covariance * a_in;
      theta_grad_in   = laplace_likelihood::theta_grad(ll_fun, theta_in, ll_args, msgs);
      obj_in     = obj_fun(a_in, theta_in);
      dir_deriv_in = grad_fun(a_in, theta_in, theta_grad_in).dot(p);
  };
  double alpha_high = alpha_init * 2;
  Eigen::VectorXd a_try(a_prev.size());
  Eigen::VectorXd theta_try(a_prev.size());
  double obj_high = 0;
  double dir_deriv_high = 0;
  update_step_vals(a_try, theta_try, theta_grad, obj_high, dir_deriv_high, alpha_high);

  // If current alpha fails, backtrack down till we find a good point
  debug::print("Begin Loop: ", 1,
    "Initial alpha: ", alpha_high,
    "g0:            ", g0.transpose().eval(),
    "theta_try:     ", theta_try.transpose().eval());
  int loop_iter = 0;
  const auto grad_tol = opt.line_search.abs_grad_threshold;
  const auto obj_tol = opt.line_search.abs_obj_threshold;
  while (alpha_high < 2) {
    // 1. Evaluate f(α) and g(α)
    update_step_vals(a_try, theta_try, theta_grad, obj_high, dir_deriv_high, alpha_high);

    debug::print("First While", 1,
          "Iter:       ", loop_iter++,
          "alpha_high: ", alpha_high,
          "obj_high:   ", obj_high,
          "deriv_high: ", dir_deriv_high,
          "deriv_init: ", dir_deriv_init,
          "theta_try:  ", theta_try.transpose());
    const bool finite_ok = std::isfinite(obj_high) && theta_try.allFinite();

    // 2. Handle numerical trouble first
    if (!finite_ok) {                      //   f or g is NaN/Inf → shrink
      alpha_high *= 0.5;
      if (alpha_high < opt.line_search.min_alpha) break;
      continue;
    }
    const bool armijo_ok = check_armijo(obj_high, obj_init,
                                        alpha_high, dir_deriv_init, opt);

    // 3. Armijo test drives all decisions
    if (armijo_ok) {                       // Armijo ✓
      const bool curve_ok  = check_wolfe_curve(dir_deriv_high, dir_deriv_init, opt);
      if (curve_ok) {
        a.swap(a_try);
        theta.swap(theta_try);
        obj_init   = obj_high;
        alpha_init = alpha_high;
        return wolfe_return::PASS;
      } else {
        alpha_high *= 2.0;
        continue;
      }
    }
    if (std::abs(dir_deriv_high) <= grad_tol              ||     // tiny slope
        std::abs(obj_high - obj_init) <= obj_tol) {               // tiny gain
      a.swap(a_try);
      theta.swap(theta_try);
      obj_init   = obj_high;
      alpha_init = alpha_high;
      return wolfe_return::GRAD_CONV;      // add a code in your enum
    }


    break;
  }
    update_step_vals(a_try, theta_try, theta_grad, obj_high, dir_deriv_high, alpha_high);

  debug::print("_______\nEnd First While: ", 1,
    "alpha_high: ", alpha_high,
    "obj_high:   ", obj_high,
    "dir_deriv_high: ", dir_deriv_high,
    "dir_deriv_init: ", dir_deriv_init,
    "theta_try:  ", theta_try.transpose());
  double alpha_success = alpha_high;
  // we *know* alpha_high is good so set that to low
  double alpha_low = 0;
  double obj_low = obj_init;
  double dir_deriv_low = dir_deriv_init;
  loop_iter = 0;
  double obj_mid = obj_low;
  double alpha_mid = alpha_low;
  double dir_deriv_mid = 0;
while (alpha_high - alpha_low > opt.line_search.min_alpha) {
  const double diff_alpha = alpha_high - alpha_low;

  alpha_mid = CubicInterp(
      alpha_low,                     /* x0  */
      -obj_low,                      /* f0  */
      -dir_deriv_low,                /* df0 */
      alpha_high,                    /* x1  */
      -obj_high,                     /* f1  */
      -dir_deriv_high,               /* df1 */
      alpha_low, alpha_high);        /* bounds */

  /* Guard against pathological cases or a NaN from the cubic.
     Fall back to bisection and keep a margin away from the ends. */
  if (!std::isfinite(alpha_mid) ||
      alpha_mid <= alpha_low || alpha_mid >= alpha_high) {
    alpha_mid = 0.5 * (alpha_low + alpha_high);
  }
  alpha_mid = std::clamp(alpha_mid, opt.line_search.min_alpha, alpha_high * 0.9);

  // --- 2. Evaluate f(α_mid) and g(α_mid) ---------------------------
  update_step_vals(a_try, theta_try, theta_grad, obj_mid, dir_deriv_mid, alpha_mid);
    debug::print("Cube: ", 1,
      "Iter:           ", loop_iter++,
      "alpha_mid:      ", alpha_mid,
      "obj_mid:        ", obj_mid,
      "dir_deriv_high: ", dir_deriv_mid,
      "alpha_low:      ", alpha_low,
      "obj_low:        ", obj_low,
      "dir_deriv_low:  ", dir_deriv_low,
      "alpha_high:     ", alpha_high,
      "obj_high:       ", obj_high,
      "dir_deriv_high: ", dir_deriv_high);
  const bool armijo_ok = check_armijo(obj_mid, obj_init,
                                      alpha_mid, dir_deriv_init, opt);
  const bool curve_ok  = check_wolfe_curve(dir_deriv_mid,
                                           dir_deriv_init,     opt);

  // --- 3. Wolfe tests & bracket update ----------------------------
  if (armijo_ok && curve_ok) {                 // (✓,✓)  accept
    a.swap(a_try);
    theta.swap(theta_try);
    obj_init   = obj_mid;
    alpha_init = alpha_mid;
    return wolfe_return::PASS;
  }

  if (!armijo_ok || obj_mid >= obj_low) {      // Armijo✗ or f↑  → shrink high
    alpha_high     = alpha_mid;
    obj_high       = obj_mid;
    dir_deriv_high = dir_deriv_mid;
    continue;
  } else {                                     // Armijo✓ & f↓
    if (dir_deriv_mid * dir_deriv_low < 0) {   // sign change → shrink high
      alpha_high     = alpha_mid;
      obj_high       = obj_mid;
      dir_deriv_high = dir_deriv_mid;
    } else {                                   // otherwise grow low
      alpha_low     = alpha_mid;
      obj_low       = obj_mid;
      dir_deriv_low = dir_deriv_mid;
    }
    continue;
  }
  if (std::abs(dir_deriv_mid) <= grad_tol              ||
      std::abs(obj_mid - obj_init) <= obj_tol) {
    a.swap(a_try);
    theta.swap(theta_try);
    obj_init   = obj_mid;
    alpha_init = alpha_mid;
    return wolfe_return::GRAD_CONV;
  }

}
  debug::print("Failed zoom: ", 1,
    "alpha_low: ", alpha_low,
    "obj_low:   ", obj_low,
    "deriv_low: ", dir_deriv_low,
    "alpha_high:", alpha_high,
    "obj_high:  ", obj_high,
    "deriv_high:", dir_deriv_high);
  // TODO: I should probably just take the largest step that armijo satisfied
  // Could not find a step that satisfied
  // accept the best good step found so far (alpha_mid)
  a.swap(a_try);
  theta.swap(theta_try);
  theta_grad   = laplace_likelihood::theta_grad(ll_fun, theta_try, ll_args, msgs);
  // We already calculated obj_mid and theta_grad so no need to recompute here
  dir_deriv_mid = grad_fun(a, theta, theta_grad).dot(p);
  const bool armijo_ok
      = check_armijo(obj_mid, obj_init, alpha_mid, dir_deriv_init, opt);
  const bool curve_ok = check_wolfe_curve(dir_deriv_mid, dir_deriv_init, opt);
  obj_init = obj_mid;
  alpha_init = alpha_mid;
  if (armijo_ok && curve_ok) {
    return wolfe_return::PASS;
  } else if (armijo_ok) {
    return wolfe_return::ARMIJO_PASS;
  } else {
    return wolfe_return::FAIL;
  }
}

/**
 * For a latent Gaussian model with hyperparameters phi and
 * latent variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approxmation. The latter is obtained by finding the mode,
 * via Newton's method, and computing the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton is a small change in
 * log marginal density. The user controls the tolerance (i.e.
 * threshold under which change is deemed small enough) and
 * maximum number of steps.
 *
 * A description of this algorithm can be found in:
 *  - (2023) Margossian, "General Adjoint-Differentiated Laplace approximation",
 *    https://arxiv.org/pdf/2306.14976.
 * Additional references include:
 *  - (2020) Margossian et al, "HMC using an adjoint-differentiated Laplace...",
 *    NeurIPS, https://arxiv.org/abs/2004.12550.
 *  - (2006) Rasmussen and Williams, "Gaussian Processes for Machine Learning",
 *    second edition, MIT Press, algorithm 3.1.
 *
 * Variables needed for the gradient or generating quantities
 * are stored by reference.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec,  InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 *
 * @return A struct containing
 * 1. lmd the log marginal density, p(y | phi)
 * 2. covariance the evaluated covariance function for the latent gaussian
 * variable
 * 3. theta a vector to store the mode
 * 4. W_r A sparse matrix containing the square root of the negative
 *    hessian, if solver 1 or 2 are used.
 * 5. L cholesky decomposition of stabilized inverse covariance
 * 6. a element in the Newton step
 * 7. l_grad the log density of the likelihood, evaluated at the mode
 *
 */
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename CovarArgs, bool InitTheta,
          require_t<is_all_arithmetic_scalar<CovarArgs>>* = nullptr>
inline auto laplace_marginal_density_est(
    LLFun&& ll_fun, LLTupleArgs&& ll_args, CovarFun&& covariance_function,
    CovarArgs&& covar_args, const laplace_options<InitTheta>& options,
    std::ostream* msgs) {
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::VectorXd;
  if constexpr (InitTheta) {
    check_nonzero_size("laplace_marginal", "initial guess", options.theta_0);
    check_finite("laplace_marginal", "initial guess", options.theta_0);
  }
  check_nonnegative("laplace_marginal", "tolerance", options.tolerance);
  check_positive("laplace_marginal", "max_num_steps", options.max_num_steps);
  check_positive("laplace_marginal", "hessian_block_size",
                 options.hessian_block_size);

  Eigen::MatrixXd covariance = stan::math::apply(
      [msgs, &covariance_function](auto&&... args) {
        return covariance_function(args..., msgs);
      },
      covar_args);
  check_square("laplace_marginal", "covariance", covariance);

  const Eigen::Index theta_size = covariance.rows();

  if (unlikely(theta_size % options.hessian_block_size != 0)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "laplace_marginal_density: The hessian size (" << theta_size
          << ", " << theta_size
          << ") is not divisible by the hessian block size ("
          << options.hessian_block_size
          << ")"
             ". Try a hessian block size such as [1, ";
      for (int i = 2; i < 12; ++i) {
        if (theta_size % i == 0) {
          msg << i << ", ";
        }
      }
      msg.str().pop_back();
      msg.str().pop_back();
      msg << "].";
      throw std::domain_error(msg.str());
    }();
  } else if (unlikely(theta_size < options.hessian_block_size)) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "laplace_marginal_density: The hessian size (" << theta_size
          << ", " << theta_size << ") is smaller than the hessian block size ("
          << options.hessian_block_size
          << "). Try a hessian block size such as [1, ";
      for (int i = 2; i < theta_size; ++i) {
        if (theta_size % i == 0) {
          msg << i << ", ";
        }
      }
      msg.str().pop_back();
      msg.str().pop_back();
      msg << "].";
      throw std::domain_error(msg.str());
    }();
  }

  auto throw_overstep = [](const auto max_num_steps) STAN_COLD_PATH {
    throw std::domain_error(
        std::string("laplace_marginal_density: max number of iterations: ")
        + std::to_string(max_num_steps) + " exceeded.");
  };
  auto ll_args_vals = value_of(ll_args);
  Eigen::VectorXd theta = [theta_size, &options]() {
    if constexpr (InitTheta) {
      return options.theta_0;
    } else {
      return Eigen::VectorXd::Zero(theta_size);
    }
  }();
  Eigen::VectorXd a_prev = Eigen::VectorXd::Zero(theta_size);
  Eigen::MatrixXd B(theta_size, theta_size);
  Eigen::VectorXd a(theta_size);
  Eigen::VectorXd b(theta_size);
  Eigen::VectorXd a_tmp(theta_size);
  Eigen::VectorXd theta_tmp(theta_size);
  // FIXME: We should use less full scope referencing here. Hard to follow
  auto obj_fun = [&](const Eigen::VectorXd& a_val, auto&& theta_val) -> double {
    return -0.5 * a_val.dot(theta_val)
           + laplace_likelihood::log_likelihood(ll_fun, theta_val, ll_args_vals,
                                                msgs);
  };
  auto grad_fun = [&](const Eigen::VectorXd& a_val, auto&& theta_val,
                      auto&& theta_grad) -> Eigen::VectorXd {
    return -covariance * a_val + covariance * theta_grad;
  };
  double objective_old = std::numeric_limits<double>::lowest();
  double objective_new = obj_fun(a_prev, theta);
  if (!std::isfinite(objective_new)) {
    throw std::domain_error(
        "laplace_marginal_density: log likelihood is not finite at initial "
        "theta and likelihood arguments.");
  }
  double step_size = 1.0;
  // NOTE: theta_grad is updated in `wolfe_line_search`
  auto theta_grad
      = laplace_likelihood::theta_grad(ll_fun, theta, ll_args, msgs);
  int iter = 0;
  if (options.solver == 1) {
    if (options.hessian_block_size == 1) {
      for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
        debug::print("======\nIter", iter++);
        auto W = laplace_likelihood::block_hessian(
            ll_fun, theta, options.hessian_block_size, ll_args, msgs);
        Eigen::VectorXd W_r(W.rows());
        Eigen::VectorXd W_vec(W.rows());
        for (Eigen::Index i = 0; i < W.rows(); i++) {
          if (W.coeff(i, i) < 0) {
            throw std::domain_error(
                "laplace_marginal_density: Hessian matrix is not positive "
                "definite");
          } else {
            W_r.coeffRef(i) = std::sqrt(W.coeff(i, i));
            W_vec.coeffRef(i) = -W.coeff(i, i);
          }
        }
        B.noalias() = MatrixXd::Identity(theta_size, theta_size)
                      + W_r.asDiagonal() * covariance * W_r.asDiagonal();
        Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
        auto L = llt_B.matrixL();
        auto LT = llt_B.matrixU();
        b.noalias() = W.diagonal().cwiseProduct(theta) + theta_grad;
        a.noalias()
            = b
              - W_r.asDiagonal()
                    * LT.solve(L.solve(W_r.cwiseProduct(covariance * b)));
      Eigen::VectorXd p = a - a_prev;
      const Eigen::VectorXd g0 = -covariance * a_prev + covariance * theta_grad;
      double g0_dir = g0.dot(p);
      double d0_dir = p.dot(W_vec.asDiagonal() * p);  // pᵀHp
      auto newton_step_size =
          (d0_dir < 0.0) ? -g0_dir / d0_dir     // concave => quadratic maximiser
                        : step_size;                 // fallback if H says ‘convex’
      debug::print("", 1,
            "Newton step size: ", newton_step_size,
            "prev step size:   ", step_size);
      step_size = newton_step_size <= 0.2 ? step_size : newton_step_size;
        objective_old = objective_new;
        auto ok = internal::wolfe_line_search(
            theta, objective_new, step_size, theta_grad, a, a_prev, ll_fun,
            obj_fun, grad_fun, covariance, ll_args, options, msgs);
        // Check for convergence or if line search failed
        debug::print("", 1,
            "Objective old: ", objective_old,
            "Objective new: ", objective_new,
            "Step size:      ", step_size);
        if (abs(objective_new - objective_old) < options.tolerance
            || (ok != wolfe_return::PASS && objective_new == objective_old)) {
          const double B_log_determinant
              = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
          // Overwrite W instead of making a new sparse matrix
          W.diagonal() = W_r;
          return laplace_density_estimates{
              objective_new - 0.5 * B_log_determinant,
              std::move(covariance),
              std::move(theta),
              std::move(W),
              Eigen::MatrixXd(L),
              std::move(a),
              std::move(theta_grad),
              Eigen::PartialPivLU<Eigen::MatrixXd>{},
              Eigen::MatrixXd(0, 0)};
        } else {
          a_prev.swap(a);
          set_zero_adjoint(ll_args);
          step_size = step_size < options.line_search.min_alpha ? 1 : step_size;
        }
      }
    } else {
      Eigen::SparseMatrix<double> W_r(theta_size, theta_size);
      Eigen::Index block_size = options.hessian_block_size;
      W_r.reserve(Eigen::VectorXi::Constant(W_r.cols(), block_size));
      const Eigen::Index n_block = W_r.cols() / block_size;
      // Prefill W_r so we only make space once
      for (Eigen::Index i = 0; i < n_block; i++) {
        for (Eigen::Index k = 0; k < block_size; k++) {
          for (Eigen::Index j = 0; j < block_size; j++) {
            W_r.insert(i * block_size + j, i * block_size + k) = 1.0;
          }
        }
      }
      W_r.makeCompressed();
      for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
        debug::print("======\nIter", iter++);
        auto W = laplace_likelihood::block_hessian(
            ll_fun, theta, options.hessian_block_size, ll_args, msgs);
        for (Eigen::Index i = 0; i < W.rows(); i++) {
          if (W.coeff(i, i) < 0) {
            throw std::domain_error(
                "laplace_marginal_density: Hessian matrix is not positive "
                "definite");
          }
        }
        block_matrix_sqrt(W_r, W, options.hessian_block_size);
        B.noalias() = MatrixXd::Identity(theta_size, theta_size)
                      + W_r * (covariance * W_r);
        Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
        if (llt_B.info() != Eigen::Success) {
          throw std::domain_error(
              "laplace_marginal_density: Cholesky failed in iteration "
              + std::to_string(i));
        }
        auto L = llt_B.matrixL();
        auto LT = llt_B.matrixU();
        b.noalias() = W * theta + theta_grad;
        a.noalias() = b - W_r * LT.solve(L.solve(W_r * (covariance * b)));
        // Simple Newton step
        objective_old = objective_new;
        Eigen::VectorXd p = a - a_prev;
        const Eigen::VectorXd g0 = -covariance * a_prev + covariance * theta_grad;
        double g0_dir = g0.dot(p);
        double d0_dir = p.dot((-W) * p);  // pᵀHp
        auto newton_step_size =
            (d0_dir < 0.0) ? -g0_dir / d0_dir     // concave => quadratic maximiser
                          : step_size;                 // fallback if H says ‘convex’
        debug::print("", 1,
            "Newton step size: ", newton_step_size,
            "prev step size:   ", step_size);
        step_size = newton_step_size <= 0.2 ? step_size : newton_step_size;

        auto ok = internal::wolfe_line_search(
            theta, objective_new, step_size, theta_grad, a, a_prev, ll_fun,
            obj_fun, grad_fun, covariance, ll_args, options, msgs);
        debug::print("", 1,
            "Objective old: ", objective_old,
            "Objective new: ", objective_new,
            "Step size:      ", step_size);
        // Check for convergence or if line search failed
        if (abs(objective_new - objective_old) < options.tolerance
            || (ok != wolfe_return::PASS && objective_new == objective_old)) {
          const double B_log_determinant
              = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
          return laplace_density_estimates{
              objective_new - 0.5 * B_log_determinant,
              std::move(covariance),
              std::move(theta),
              std::move(W_r),
              Eigen::MatrixXd(L),
              std::move(a),
              std::move(theta_grad),
              Eigen::PartialPivLU<Eigen::MatrixXd>{},
              Eigen::MatrixXd(0, 0)};
        } else {
          a_prev.swap(a);
          set_zero_adjoint(ll_args);
        }
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 2) {
    Eigen::MatrixXd K_root
        = covariance.template selfadjointView<Eigen::Lower>().llt().matrixL();
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      debug::print("======\nIter", iter++);
      auto W = laplace_likelihood::block_hessian(
          ll_fun, theta, options.hessian_block_size, ll_args, msgs);
      B.noalias() = MatrixXd::Identity(theta_size, theta_size)
                    + K_root.transpose() * W * K_root;
      Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt_B(B);
      if (llt_B.info() != Eigen::Success) {
        throw std::domain_error(
            "laplace_marginal_density: Cholesky failed in iteration "
            + std::to_string(i));
      }
      auto L = llt_B.matrixL();
      auto LT = llt_B.matrixU();
      b.noalias() = W * theta + theta_grad;
      a.noalias()
          = K_root.transpose().template triangularView<Eigen::Upper>().solve(
              LT.solve(L.solve(K_root.transpose() * b)));
      objective_old = objective_new;
      Eigen::VectorXd p = a - a_prev;
      const Eigen::VectorXd g0 = -covariance * a_prev + covariance * theta_grad;
      double g0_dir = g0.dot(p);
      double d0_dir = p.dot((-W) * p);  // pᵀHp
      auto newton_step_size =
          (d0_dir < 0.0) ? -g0_dir / d0_dir     // concave => quadratic maximiser
                        : step_size;                 // fallback if H says ‘convex’
      debug::print("", 1,
            "Newton step size: ", newton_step_size,
            "prev step size:   ", step_size);
      step_size = newton_step_size <= 0.2 ? step_size : newton_step_size;

      auto ok = internal::wolfe_line_search(
          theta, objective_new, step_size, theta_grad, a, a_prev, ll_fun,
          obj_fun, grad_fun, covariance, ll_args, options, msgs);
      // Check for convergence or if line search failed
        debug::print("", 1,
            "Objective old: ", objective_old,
            "Objective new: ", objective_new,
            "Step size:      ", step_size);
      if (abs(objective_new - objective_old) < options.tolerance
          || (ok != wolfe_return::PASS && objective_new == objective_old)) {
        const double B_log_determinant
            = 2.0 * llt_B.matrixLLT().diagonal().array().log().sum();
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W),
            std::move(Eigen::MatrixXd(L)),
            std::move(a),
            std::move(theta_grad),
            Eigen::PartialPivLU<Eigen::MatrixXd>{},
            std::move(K_root)};
      } else {
        a_prev = a;
        set_zero_adjoint(ll_args);
      }
    }
    throw_overstep(options.max_num_steps);
  } else if (options.solver == 3) {
    for (Eigen::Index i = 0; i <= options.max_num_steps; i++) {
      debug::print("======\nIter", iter++);
      auto W = laplace_likelihood::block_hessian(
          ll_fun, theta, options.hessian_block_size, ll_args, msgs);
      Eigen::PartialPivLU<Eigen::MatrixXd> LU(
          MatrixXd::Identity(theta_size, theta_size) + covariance * W);
      // L on lower and U on upper triangular
      b.noalias() = W * theta + theta_grad;
      a.noalias() = b - W * LU.solve(covariance * b);
      objective_old = objective_new;
        Eigen::VectorXd p = a - a_prev;
        const Eigen::VectorXd g0 = -covariance * a_prev + covariance * theta_grad;
        double g0_dir = g0.dot(p);
        double d0_dir = p.dot((-W) * p);  // pᵀHp
        auto newton_step_size =
            (d0_dir < 0.0) ? -g0_dir / d0_dir     // concave => quadratic maximiser
                          : step_size;                 // fallback if H says ‘convex’
      debug::print("", 1,
            "Newton step size: ", newton_step_size,
            "prev step size:   ", step_size);
      step_size = newton_step_size <= 0.2 ? step_size : newton_step_size;
      auto ok = internal::wolfe_line_search(
          theta, objective_new, step_size, theta_grad, a, a_prev, ll_fun,
          obj_fun, grad_fun, covariance, ll_args, options, msgs);
        debug::print("", 1,
            "Objective old: ", objective_old,
            "Objective new: ", objective_new,
            "Step size:      ", step_size);
      // Check for convergence or if line search failed
      if (abs(objective_new - objective_old) < options.tolerance
          || (ok != wolfe_return::PASS && objective_new == objective_old)) {
        // TODO(Charles): There has to be a simple trick for this
        const double B_log_determinant = log(LU.determinant());
        return laplace_density_estimates{
            objective_new - 0.5 * B_log_determinant,
            std::move(covariance),
            std::move(theta),
            std::move(W),
            Eigen::MatrixXd(0, 0),
            std::move(a),
            std::move(theta_grad),
            std::move(LU),
            Eigen::MatrixXd(0, 0)};
      } else {
        a_prev = a;
        set_zero_adjoint(ll_args);
        step_size = step_size < 1e-3 ? 1 : step_size;
      }
    }
    throw_overstep(options.max_num_steps);
  }
  throw std::domain_error(
      std::string("You chose a solver (") + std::to_string(options.solver)
      + ") that is not valid. Please choose either 1, 2, or 3.");
}
}  // namespace internal
/**
 * For a latent Gaussian model with global parameters phi, latent
 * variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approxmation. The latter is obtained by finding the mode,
 * using a custom Newton method, and the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton is a small change in
 * log marginal density. The user controls the tolerance (i.e.
 * threshold under which change is deemed small enough) and
 * maximum number of steps.
 *
 * Wrapper for when the hyperparameters are passed as a double.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec, InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 * @return the log maginal density, p(y | phi)
 */
template <
    typename LLFun, typename LLTupleArgs, typename CovarFun, typename CovarArgs,
    bool InitTheta,
    require_t<is_all_arithmetic_scalar<CovarArgs, LLTupleArgs>>* = nullptr>
inline double laplace_marginal_density(
    LLFun&& ll_fun, LLTupleArgs&& ll_args, CovarFun&& covariance_function,
    CovarArgs&& covar_args, const laplace_options<InitTheta>& options,
    std::ostream* msgs) {
  return internal::laplace_marginal_density_est(
             std::forward<LLFun>(ll_fun), std::forward<LLTupleArgs>(ll_args),
             std::forward<CovarFun>(covariance_function),
             std::forward<CovarArgs>(covar_args), options, msgs)
      .lmd;
}

namespace internal {

/**
 * Collects the adjoints from the input and adds them to the output.
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `arithmetic` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_all_arithmetic_scalar<Input>>* = nullptr>
inline void collect_adjoints(Output&& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<double, -1, 1>> input_map(input_i.data(),
                                                             input_i.size());
          output_map.array() += input_map.array();
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += input_i.array();
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += input_i;
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}
/**
 * Base case for zero sized tuples
 */
template <bool ZeroInput = false>
inline void constexpr copy_compute_s2(const std::tuple<>& output,
                                      const std::tuple<>& input) noexcept {}

/**
 * Copies the adjoints from the input to the output, scaling them by 0.5.
 * @tparam ZeroInput If true, the adjoints of the input will be set to zero
 * @tparam Output A tuple or type where all scalar types are `arithmetic` types
 * @tparam Input A tuple or type where all scalar types are `var` types
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <bool ZeroInput = false, typename Output, typename Input,
          require_t<is_all_arithmetic_scalar<Output>>* = nullptr,
          require_t<is_any_var_scalar<Input>>* = nullptr>
inline void copy_compute_s2(Output&& output, Input&& input) {
  return iter_tuple_nested(
      [](auto&& output_i, auto&& input_i) {
        using output_i_t = std::decay_t<decltype(output_i)>;
        if constexpr (is_std_vector_v<output_i_t>) {
          Eigen::Map<Eigen::Matrix<double, -1, 1>> output_map(output_i.data(),
                                                              output_i.size());
          Eigen::Map<Eigen::Matrix<var, -1, 1>> input_map(input_i.data(),
                                                          input_i.size());
          output_map.array() += 0.5 * input_map.adj().array();
          if constexpr (ZeroInput) {
            input_map.adj().setZero();
          }
        } else if constexpr (is_eigen_v<output_i_t>) {
          output_i.array() += 0.5 * input_i.adj().array();
          if constexpr (ZeroInput) {
            input_i.adj().setZero();
          }
        } else if constexpr (is_stan_scalar_v<output_i_t>) {
          output_i += (0.5 * input_i.adj());
          if constexpr (ZeroInput) {
            input_i.adj() = 0;
          }
        } else {
          static_assert(
              sizeof(std::decay_t<output_i_t>*) == 0,
              "INTERNAL ERROR:(laplace_marginal_lpdf) set_zero_adjoints was "
              "not able to deduce the actions needed for the given type. "
              "This is an internal error, please report it: "
              "https://github.com/stan-dev/math/issues");
        }
      },
      std::forward<Output>(output), std::forward<Input>(input));
}

template <typename T>
inline constexpr decltype(auto) filter_var_scalar_types(T&& t) {
  return stan::math::filter_map<is_any_var_scalar>(
      [](auto&& arg) -> decltype(auto) {
        return std::forward<decltype(arg)>(arg);
      },
      std::forward<T>(t));
}
/**
 * Creates an arena type from the input with initialized with zeros
 * @tparam Input Possibly a tuple, std::vector, Eigen type, or scalar
 * @param input The input to be converted to an arena type
 */
template <typename Input>
inline constexpr auto make_zeroed_arena(Input&& input) {
  if constexpr (is_tuple_v<Input>) {
    return stan::math::filter_map<is_any_var_scalar>(
        [](auto&& output_i) { return make_zeroed_arena(output_i); }, input);
  } else if constexpr (is_std_vector_v<Input>) {
    if constexpr (!is_var_v<value_type_t<Input>>) {
      const auto output_size = input.size();
      arena_t<std::vector<decltype(make_zeroed_arena(input[0]))>> ret;
      ret.reserve(output_size);
      for (Eigen::Index i = 0; i < output_size; ++i) {
        ret.push_back(make_zeroed_arena(input[i]));
      }
      return ret;
    } else {
      return arena_t<std::vector<double>>(input.size(), 0.0);
    }
  } else if constexpr (is_eigen_v<Input>) {
    return arena_t<promote_scalar_t<double, Input>>(
        plain_type_t<promote_scalar_t<double, Input>>::Zero(input.rows(),
                                                            input.cols()));
  } else if constexpr (is_var<Input>::value) {
    return static_cast<double>(0.0);
  }
}

/**
 * Used in reverse pass to collect adjoints to the output
 * @tparam Output A tuple or type where all scalar types are `var` types
 * @tparam Input A tuple or type where all scalar types are `arithmetic` types
 * @param output The output to which the adjoints will be added
 * @param ret The vari object containing the adjoint to be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input>
inline void collect_adjoints(Output&& output, const vari* ret, Input&& input) {
  if constexpr (is_tuple_v<Output>) {
    static_assert(sizeof(std::decay_t<Output>*) == 0,
                  "INTERNAL ERROR:(laplace_marginal_lpdf) "
                  "Accumulate Adjoints called on a tuple, but tuples cannot be "
                  "on the reverse mode stack! "
                  "This is an internal error, please report it: "
                  "https://github.com/stan-dev/math/issues");
  } else if constexpr (is_std_vector_v<Output>) {
    if constexpr (!is_var_v<value_type_t<Output>>) {
      const auto output_size = output.size();
      for (std::size_t i = 0; i < output_size; ++i) {
        collect_adjoints(output[i], ret, input[i]);
      }
    } else {
      Eigen::Map<Eigen::Matrix<var, -1, 1>> output_map(output.data(),
                                                       output.size());
      Eigen::Map<const Eigen::Matrix<double, -1, 1>> input_map(input.data(),
                                                               input.size());
      output_map.array().adj() += ret->adj_ * input_map.array();
    }
  } else if constexpr (is_eigen_v<Output>) {
    output.adj().array() += ret->adj_ * input.array();
  } else if constexpr (is_var_v<Output>) {
    output.adj() += ret->adj_ * input;
  }
}

/**
 * Collects adjoints from a tuple or std::vector of tuples
 * @tparam Output A tuple or std::vector of tuples where all scalar types are
 * `var` types
 * @tparam Input A tuple or std::vector of tuples where all scalar types are
 * `arithmetic` types
 * @param ret The vari object containing the adjoint to be added
 * @param output The output to which the adjoints will be added
 * @param input The input from which the adjoints will be collected
 */
template <typename Output, typename Input>
inline void reverse_pass_collect_adjoints(var ret, Output&& output,
                                          Input&& input) {
  if constexpr (is_tuple_v<Output>) {
    stan::math::for_each(
        [ret](auto&& inner_arg, auto&& inner_input) mutable {
          reverse_pass_collect_adjoints(
              ret, std::forward<decltype(inner_arg)>(inner_arg),
              std::forward<decltype(inner_input)>(inner_input));
        },
        std::forward<Output>(output), std::forward<Input>(input));
  } else if constexpr (is_std_vector_containing_tuple_v<Output>) {
    for (std::size_t i = 0; i < output.size(); ++i) {
      reverse_pass_collect_adjoints(ret, output[i], input[i]);
    }
  } else {
    reverse_pass_callback(
        [vi = ret.vi_, arg_arena = to_arena(std::forward<Output>(output)),
         input_arena = to_arena(std::forward<Input>(input))]() mutable {
          collect_adjoints(arg_arena, vi, input_arena);
        });
  }
}
}  // namespace internal
/**
 * For a latent Gaussian model with global parameters phi, latent
 * variables theta, and observations y, this function computes
 * an approximation of the log marginal density, p(y | phi).
 * This is done by marginalizing out theta, using a Laplace
 * approxmation. The latter is obtained by finding the mode,
 * using a custom Newton method, and the Hessian of the likelihood.
 *
 * The convergence criterion for the Newton is a small change in
 * the log marginal density. The user controls the tolerance (i.e.
 * threshold under which change is deemed small enough) and
 * maximum number of steps.
 *
 * Wrapper for when the global parameter is passed as a double.
 *
 * @tparam LLFun Type with a valid `operator(ThetaVec,  InnerLLTupleArgs)`
 * where `InnerLLTupleArgs` are the elements of `LLTupleArgs`
 * @tparam LLTupleArgs A tuple whose elements follow the types required for
 * `LLFun`
 * \laplace_common_template_args
 * @param[in] ll_fun A log likelihood functor
 * @param[in] ll_args Tuple containing parameters for `LLFun`
 * \laplace_common_args
 * @param[in] options A set of options for tuning the solver
 * \msg_arg
 * @return the log maginal density, p(y | phi)
 */
template <typename LLFun, typename LLTupleArgs, typename CovarFun,
          typename CovarArgs, bool InitTheta,
          require_t<is_any_var_scalar<LLTupleArgs, CovarArgs>>* = nullptr>
inline auto laplace_marginal_density(const LLFun& ll_fun, LLTupleArgs&& ll_args,
                                     CovarFun&& covariance_function,
                                     CovarArgs&& covar_args,
                                     const laplace_options<InitTheta>& options,
                                     std::ostream* msgs) {
  auto covar_args_refs = to_ref(std::forward<CovarArgs>(covar_args));
  auto ll_args_refs = to_ref(std::forward<LLTupleArgs>(ll_args));
  // Solver 1, 2, 3
  constexpr bool ll_args_contain_var = is_any_var_scalar<LLTupleArgs>::value;
  auto partial_parm = internal::make_zeroed_arena(ll_args_refs);
  auto covar_args_adj = internal::make_zeroed_arena(covar_args_refs);
  double lmd = 0.0;
  {
    nested_rev_autodiff nested;

    // Make one hard copy here
    using laplace_likelihood::internal::conditional_copy_and_promote;
    using laplace_likelihood::internal::COPY_TYPE;
    auto ll_args_copy
        = conditional_copy_and_promote<is_any_var_scalar, var, COPY_TYPE::DEEP>(
            ll_args_refs);

    auto md_est = internal::laplace_marginal_density_est(
        ll_fun, ll_args_copy, covariance_function, value_of(covar_args_refs),
        options, msgs);

    // Solver 1, 2
    arena_t<Eigen::MatrixXd> R(md_est.theta.size(), md_est.theta.size());
    // Solver 3
    arena_t<Eigen::MatrixXd> LU_solve_covariance;
    // Solver 1, 2, 3
    arena_t<Eigen::VectorXd> s2(md_est.theta.size());

    // Return references to var types
    auto ll_args_filter = internal::filter_var_scalar_types(ll_args_copy);
    stan::math::for_each(
        [](auto&& output_i, auto&& ll_arg_i) {
          if (is_any_var_scalar_v<decltype(ll_arg_i)>) {
            internal::collect_adjoints<true>(output_i, ll_arg_i);
          }
        },
        partial_parm, ll_args_filter);
    if (options.solver == 1) {
      if (options.hessian_block_size == 1) {
        // TODO(Steve): Solve without casting from sparse to dense
        Eigen::MatrixXd tmp
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r.toDense());
        R = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * md_est.covariance);
        if constexpr (!ll_args_contain_var) {
          s2.deep_copy(
              (0.5
               * (md_est.covariance.diagonal() - (C.transpose() * C).diagonal())
                     .cwiseProduct(laplace_likelihood::third_diff(
                         ll_fun, md_est.theta, value_of(ll_args_copy), msgs))));
        } else {
          arena_t<Eigen::MatrixXd> A = md_est.covariance - C.transpose() * C;
          auto s2_tmp = laplace_likelihood::compute_s2(
              ll_fun, md_est.theta, A, options.hessian_block_size, ll_args_copy,
              msgs);
          s2.deep_copy(s2_tmp);
          internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
        }

      } else {
        Eigen::MatrixXd tmp
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r.toDense());
        R = tmp.transpose() * tmp;
        arena_t<Eigen::MatrixXd> C
            = md_est.L.template triangularView<Eigen::Lower>().solve(
                md_est.W_r * md_est.covariance);
        arena_t<Eigen::MatrixXd> A = md_est.covariance - C.transpose() * C;
        auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                     options.hessian_block_size,
                                                     ll_args_copy, msgs);
        s2.deep_copy(s2_tmp);
        internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
      }
    } else if (options.solver == 2) {
      R = md_est.W_r
          - md_est.W_r * md_est.K_root
                * md_est.L.transpose()
                      .template triangularView<Eigen::Upper>()
                      .solve(
                          md_est.L.template triangularView<Eigen::Lower>()
                              .solve(md_est.K_root.transpose() * md_est.W_r));

      arena_t<Eigen::MatrixXd> C
          = md_est.L.template triangularView<Eigen::Lower>().solve(
              md_est.K_root.transpose());
      auto s2_tmp = laplace_likelihood::compute_s2(
          ll_fun, md_est.theta, (C.transpose() * C).eval(),
          options.hessian_block_size, ll_args_copy, msgs);
      s2.deep_copy(s2_tmp);
      internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
    } else {  // options.solver with LU decomposition
      LU_solve_covariance = md_est.LU.solve(md_est.covariance);
      R = md_est.W_r - md_est.W_r * LU_solve_covariance * md_est.W_r;
      arena_t<Eigen::MatrixXd> A
          = md_est.covariance
            - md_est.covariance * md_est.W_r * LU_solve_covariance;
      auto s2_tmp = laplace_likelihood::compute_s2(ll_fun, md_est.theta, A,
                                                   options.hessian_block_size,
                                                   ll_args_copy, msgs);
      s2.deep_copy(s2_tmp);
      internal::copy_compute_s2<true>(partial_parm, ll_args_filter);
    }
    lmd = md_est.lmd;
    if constexpr (is_any_var_scalar_v<scalar_type_t<CovarArgs>>) {
      [&covar_args_refs, &covar_args_adj, &md_est, &R, &s2,
       &covariance_function, &msgs]() mutable {
        const nested_rev_autodiff nested;
        auto covar_args_copy
            = laplace_likelihood::internal::conditional_copy_and_promote<
                is_any_var_scalar, var,
                laplace_likelihood::internal::COPY_TYPE::DEEP>(covar_args_refs);

        var_value<Eigen::MatrixXd> K_var = to_var_value(stan::math::apply(
            [&covariance_function, &msgs](auto&&... args) {
              return covariance_function(args..., msgs);
            },
            covar_args_copy));
        arena_t<Eigen::MatrixXd> K_adj_arena
            = 0.5 * md_est.a * md_est.a.transpose() - 0.5 * R
              + s2 * md_est.theta_grad.transpose()
              - (R * (K_var.val() * s2)) * md_est.theta_grad.transpose();
        var Z = make_callback_var(0.0, [K_var, K_adj_arena](auto&& vi) mutable {
          K_var.adj().array() += vi.adj() * K_adj_arena.array();
        });
        grad(Z.vi_);
        auto covar_args_filter
            = internal::filter_var_scalar_types(covar_args_copy);
        internal::collect_adjoints(covar_args_adj, covar_args_filter);
      }();
    }
    if constexpr (ll_args_contain_var) {
      arena_t<Eigen::VectorXd> v;
      if (options.solver == 1 || options.solver == 2) {
        v = md_est.covariance * s2
            - md_est.covariance * R * md_est.covariance * s2;
      } else {
        v = LU_solve_covariance * s2;
      }
      laplace_likelihood::diff_eta_implicit(ll_fun, v, md_est.theta,
                                            ll_args_copy, msgs);
      internal::collect_adjoints<true>(partial_parm, ll_args_filter);
    }
  }
  var ret(lmd);
  if constexpr (is_any_var_scalar_v<CovarArgs>) {
    auto covar_args_filter = internal::filter_var_scalar_types(covar_args_refs);
    internal::reverse_pass_collect_adjoints(ret, covar_args_filter,
                                            covar_args_adj);
  }
  if constexpr (ll_args_contain_var) {
    auto ll_args_filter = internal::filter_var_scalar_types(ll_args_refs);
    internal::reverse_pass_collect_adjoints(ret, ll_args_filter, partial_parm);
  }
  return ret;
}

}  // namespace math
}  // namespace stan

#endif
