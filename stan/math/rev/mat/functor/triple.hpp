#ifndef STAN_MATH_PRIM_MAT_FUN_triple_HPP
#define STAN_MATH_PRIM_MAT_FUN_triple_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <stan/math/rev/core.hpp>
#include <iostream>  // TEST

/**
 * TEST: this function is here for testing purposes. 
 * Notably, it demonstrates using the Jacobian function
 * inside the chain() method causes the code to crash
 * when calling the grad() method.
 */

namespace stan {
  namespace math {

    struct triple_functor {
    public:
      triple_functor() { };

      /**
       * Multiply matrix by three.
       */  
      template <typename T>
      inline
      Eigen::Matrix<T, Eigen::Dynamic, 1>
      operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1> x) const {
        return 3.0 * x;
      }
    };

    namespace {
      template <typename T>
      class triple_vari_alloc : public chainable_alloc {
      private:
        inline void compute(const Eigen::VectorXd& x) {
          triple_functor f;
          Eigen::VectorXd y_dbl = f(x);

          for (int i = 0; i < y_dbl.size(); i++)
            y_(i) = var(new vari(y_dbl(i), false));
      }

      public:
        triple_vari_alloc(const Eigen::Matrix<T,
                            Eigen::Dynamic, 1>& x)
          : x_(x), y_(3) {
          compute(value_of(x));
         }

         Eigen::Matrix<T, Eigen::Dynamic, 1> x_;
         Eigen::Matrix<T, Eigen::Dynamic, 1> y_;
      };

      template <typename T>
      class triple_vari : public vari {
      protected:
        inline void chain(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
                          const Eigen::VectorXd adjY) {
         // Implementation of Jacobian for testing purposes:
         // the Jacobian does not contribute to the result,
         // but calling it (eventually) causes the code to 
         // break.
         Eigen::MatrixXd J;
         Eigen::VectorXd fvec;
         triple_functor f;
         stan::math::jacobian(f, value_of(x), fvec, J);

         // Compute the adjoints (these three lines are
         // the only ones that should appear in the chain
         // method.
         x(0).vi_->adj_ += 3 * adjY(0);
         x(1).vi_->adj_ += 3 * adjY(1);         
         x(2).vi_->adj_ += 3 * adjY(2);
        }

      public:
        triple_vari(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
          : vari(0.0) {
          impl_ = new triple_vari_alloc<T>(x);
        }

        virtual void chain() {
          Eigen::VectorXd adjY(impl_->y_.size());
          for (int i = 0; i < adjY.rows(); i++)
            adjY(i) = impl_->y_(i).vi_->adj_;

          chain(impl_->x_, adjY);
        }

        triple_vari_alloc<T> *impl_;
      };

      /**
       * Multiplies a matrix by three.
       */
      template <typename T>
      inline Eigen::Matrix<T, Eigen::Dynamic, 1>
      triple(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
        triple_vari<T> *baseVari
          = new triple_vari<T>(x);

        return baseVari->impl_->y_;
      }
    }  
  }   
}

#endif
