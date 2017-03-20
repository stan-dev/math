#ifndef STAN_MATH_REV_MAT_FUNCTOR_JACOBIAN_TEST_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_JACOBIAN_TEST_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stdexcept>
#include <vector>
#include <iostream>  // TEST

namespace stan {
  namespace math {

    template <typename F>
    void
    jacobian_test(const F& f,
                  const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
                  Eigen::Matrix<double, Eigen::Dynamic, 1>& fx,
                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& J) {
      std::cout << "Running Jacobian_test ..." << std::endl;  // TEST
      using Eigen::Matrix;
      using Eigen::Dynamic;
      start_nested();
      try {
        Matrix<var, Dynamic, 1> x_var(x.size());
        for (int k = 0; k < x.size(); ++k)
          x_var(k) = x(k);

        // Matrix<var, Dynamic, 1> fx_var = f(x_var);  // BUG: causes error
        /*
        std::cout << "fx_var: ";
        for (int i = 0; i < fx_var.size(); i++)
          std::cout << fx_var(i).val() << ":" << fx_var(i).adj() << " ";
        std::cout << std::endl; */

        // TEST: This alternative also causes the bug (lines 36 - 43)
        /*
        Eigen::MatrixXd fx_dbl = f(x);
        Matrix<var, Dynamic, 1> fx_var(fx_dbl.size());
        fx_var(0) = fx_dbl(0);
        fx_var(1) = fx_dbl(1);  // this line causes an error
        // TEST: this loop causes an error
        for (int k = 0; k < fx_dbl.size(); ++k)
          fx_var(k) = fx_dbl(k); */

        // fx.resize(fx_var.size());
        // for (int i = 0; i < fx_var.size(); ++i)
          // fx(i) = fx_var(i).val();
        // J.resize(fx_var.size(), x.size());
        /*
        for (int i = 0; i < fx_var.size(); ++i) {
           if (i > 0)
             set_zero_all_adjoints_nested();
          grad(fx_var(i).vi_);
          for (int k = 0; k < x.size(); ++k)
               J(i, k) = x_var(k).adj();
        } */
      } catch (const std::exception& e) {
        recover_memory_nested();
        throw;
      }
      recover_memory_nested();
      J.resize(2, 3);  // TEST
    }

  }
}
#endif
