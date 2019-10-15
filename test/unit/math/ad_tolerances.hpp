#ifndef TEST_UNIT_MATH_AD_TOLERANCES_HPP
#define TEST_UNIT_MATH_AD_TOLERANCES_HPP

namespace stan {
namespace test {

/**
 * Simple struct to hold the complete set of tolerances used to test a
 * function.  The default constructor uses default values for all
 * tolerances.  These are 1e-8 for values, 1e-4 for first derivatives,
 * 1e-3 for second derivatives, and 1e-2 for third derivatives.  The
 * names begin with the functional being evaluated, include an `fvar`
 * if the function is implemented using only forward mode, and end
 * with the quantity being calculated; for example,
 * `hessian_fvar_grad_` is the gradient calculated by the `hessian`
 * function using forward-mode autodiff.
 *
 * `gradient_val_`: 1e-8;  `gradient_grad_`: 1e-4
 *
 * `gradient_fvar_val_`: 1e-8;  `gradient_fvar_grad_`: 1e-4
 *
 * `hessian_val_` : 1e-8; `hessian_grad_`: 1e-4; `hessian_hessian_`: 1e-3
 *
 * `hessian_fvar_val_` : 1e-8; `hessian_fvar_grad_`: 1e-4;
 * `hessian_fvar_hessian_`: 1e-3
 *
 * `grad_hessian_val_` : 1e-8; `grad_hessian_hessian_`: 1e-3;
 * `grad_hessian_grad_hessian_`: 1e-2
 */
struct ad_tolerances {
  double gradient_val_;
  double gradient_grad_;
  double gradient_fvar_val_;
  double gradient_fvar_grad_;
  double hessian_val_;
  double hessian_grad_;
  double hessian_hessian_;
  double hessian_fvar_val_;
  double hessian_fvar_grad_;
  double hessian_fvar_hessian_;
  double grad_hessian_val_;
  double grad_hessian_hessian_;
  double grad_hessian_grad_hessian_;
  ad_tolerances()
      : gradient_val_(1e-8),
        gradient_grad_(1e-4),

        gradient_fvar_val_(1e-8),
        gradient_fvar_grad_(1e-4),

        hessian_val_(1e-8),
        hessian_grad_(1e-4),
        hessian_hessian_(1e-3),

        hessian_fvar_val_(1e-8),
        hessian_fvar_grad_(1e-4),
        hessian_fvar_hessian_(1e-3),

        grad_hessian_val_(1e-8),
        grad_hessian_hessian_(1e-3),
        grad_hessian_grad_hessian_(1e-2) {}
};

}  // namespace test
}  // namespace stan
#endif
