#ifndef TEST_UNIT_MATH_MIX_MAT_UTIL_AUTODIFF_TESTER_HPP
#define TEST_UNIT_MATH_MIX_MAT_UTIL_AUTODIFF_TESTER_HPP

#include <stan/math/mix/mat.hpp>
#include <test/unit/math/mix/mat/seq_reader.hpp>
#include <test/unit/math/mix/mat/seq_writer.hpp>
#include <gtest/gtest.h>
#include <stdexcept>

namespace stan {
  namespace math {
    namespace test {

      // For an example of how to use this tester, see
      // test/unit/math/mix/core/operator_addition_test.cpp

      /**
       * Test that scalars x1 and x2 are within the specified
       * tolerance, with identity behavior for infinite and NaN
       * values.
       */
      template <typename T1, typename T2>
      void expect_near(const T1& x1, const T2& x2,
                       double tol = 1e-9) {
        if (is_nan(x1) || is_nan(x2))
          EXPECT_TRUE(is_nan(x1) && is_nan(x2));
        else if (is_inf(x1) || is_inf(x2))
          EXPECT_EQ(x1, x2);
        else
          EXPECT_NEAR(x1, x2, tol);
      }


      /**
       * Tests that matrices (or vectors) x1 and x2 are same size and
       * have near values up to specified tolerance.
       */
      template <typename T, int R, int C>
      void expect_near(const Eigen::Matrix<T, R, C>& x1,
                       const Eigen::Matrix<T, R, C>& x2,
                       double tol = 1e-7) {
        EXPECT_EQ(x1.rows(), x2.rows());
        EXPECT_EQ(x1.cols(), x2.cols());
        for (int i = 0; i < x1.size(); ++i)
          expect_near(x1(i), x2(i), tol);
      }


      /**
       * Tests that the function f applied to the argument x yields
       * the expected value fx (with scalars as double).
       */
      template <typename F>
      void test_value(const F& f, const Eigen::VectorXd& x, double fx) {
        if (is_nan(fx))
          EXPECT_TRUE(is_nan(f(x)));
        else
          expect_near(fx, f(x));
      }


      /**
       * Tests that the function f applied to the argument x yields
       * the value fx and the correct first-order derivatives as
       * calaculated with the gradient functional using var.
       */
      template <typename F>
      void test_gradient(const F& f, const Eigen::VectorXd& x, double fx) {
        Eigen::VectorXd grad_ad;
        double fx_ad;
        gradient<F>(f, x, fx_ad, grad_ad);
        expect_near(fx, fx_ad);
        Eigen::VectorXd grad_fd;
        double fx_fd;
        finite_diff_gradient(f, x, fx_fd, grad_fd);
        expect_near(grad_fd, grad_ad);
      }


      /**
       * Tests that the function f applied to the argument x yields
       * the value fx and correct first-order derivatives as
       * calculated by the gradient functionional using fvar<double>
       * scalars.
       */
      template <typename F>
      void test_gradient_fvar(const F& f, const Eigen::VectorXd& x, double fx) {
        Eigen::VectorXd grad_ad;
        double fx_ad;
        gradient<double, F>(f, x, fx_ad, grad_ad);
        expect_near(fx, fx_ad);
        Eigen::VectorXd grad_fd;
        double fx_fd;
        finite_diff_gradient(f, x, fx_fd, grad_fd);
        expect_near(grad_fd, grad_ad);
      }


      /**
       * Tests that the function f applied to the argument x yields
       * the value fx and correct first- and second-order derivatives
       * as calculated by the hessian functional using fvar<var>
       * scalars.
       */
      template <typename F>
      void test_hessian_fvar(const F& f, const Eigen::VectorXd& x, double fx) {
        double fx_ad;
        Eigen::VectorXd grad_ad;
        Eigen::MatrixXd H_ad;
        hessian<double, F>(f, x, fx_ad, grad_ad, H_ad);
        expect_near(fx, fx_ad);
        double fx_fd;
        Eigen::VectorXd grad_fd;
        Eigen::MatrixXd H_fd;
        finite_diff_hessian(f, x, fx_fd, grad_fd, H_fd);
        expect_near(grad_fd, grad_ad);
        expect_near(H_fd, H_ad);
      }


      /**
       * Tests that the function f applied to the argument x yields
       * the value fx and correct first- and second-order derivatives
       * as calculated by the hessian functional using
       * fvar<fvar<double>> scalars.
       */
      template <typename F>
      void test_hessian(const F& f, const Eigen::VectorXd& x, double fx) {
        double fx_ad;
        Eigen::VectorXd grad_ad;
        Eigen::MatrixXd H_ad;
        hessian<F>(f, x, fx_ad, grad_ad, H_ad);
        expect_near(fx, fx_ad);
        double fx_fd;
        Eigen::VectorXd grad_fd;
        Eigen::MatrixXd H_fd;
        finite_diff_hessian(f, x, fx_fd, grad_fd, H_fd);
        expect_near(grad_fd, grad_ad);
        expect_near(H_fd, H_ad);
      }


      /**
       * Tests that the function f applied to the argument x yields
       * the value fx and correct first-, second-, and third-order
       * derivatives as calculated by the hessian functional using
       * fvar<fvar<var>> scalars.
       */
      template <typename F>
      void test_grad_hessian(const F& f, const Eigen::VectorXd& x, double fx) {
        double fx_ad;
        Eigen::MatrixXd H_ad;
        std::vector<Eigen::MatrixXd> grad_H_ad;
        grad_hessian(f, x, fx_ad, H_ad, grad_H_ad);
        expect_near(fx, fx_ad);
        double fx_fd;
        Eigen::MatrixXd H_fd;
        std::vector<Eigen::MatrixXd> grad_H_fd;
        grad_hessian(f, x, fx_fd, H_fd, grad_H_fd);
        expect_near(H_fd, H_ad);
        EXPECT_EQ(x.size(), grad_H_fd.size());
        for (size_t i = 0; i < grad_H_fd.size(); ++i)
          expect_near(grad_H_fd[i], grad_H_ad[i]);
      }


      // test value and derivative in all functionals vs. finite diffs
      template <typename F>
      void test_functor(const F& f, const Eigen::VectorXd& x, double fx) {
        test_value(f, x, fx);

        // finite diffs can't handle infinity
        if (is_inf(fx)) return;
        for (int i = 0; i < x.size(); ++i)
          if (is_inf(x(i)))
            return;

        test_gradient(f, x, fx);
        test_gradient_fvar(f, x, fx);
        test_hessian(f, x, fx);
        test_hessian_fvar(f, x, fx);
        // test_grad_hessian(f, x, fx);
      }


      /**
       * Structure to adapt a two-argument function specified as a
       * class with a static apply(,) method that returns a scalar to
       * a function that operates on an Eigen vector with templated
       * scalar type and returns a scalar of the same type.
       *
       * <p>It works by adapting the two-argument function to be a
       * vector function with zero, one or both arguments being
       * instantiated to doubles and the remaining arguments being
       * templated on the functor argument scalar type.
       *
       * @tparam F class with static apply(,) method
       * @tparam T1 type of first argument with double scalars
       * @tparam T2 type of second argument with double scalars
       */
      template <typename F, typename T1, typename T2>
      struct binder_binary {
        T1 x1_;
        T2 x2_;
        bool fixed1_;
        bool fixed2_;

        binder_binary(const T1& x1, const T2& x2)
          : x1_(x1), x2_(x2), fixed1_(false), fixed2_(false) { }

        template <typename T>
        T operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
          if (fixed1_ && fixed2_) {
            return F::apply(x1_, x2_);
          } else if (fixed1_ && !fixed2_) {
            seq_reader<T> r(theta);
            return F::apply(x1_, r.read(x2_));
          } else if (!fixed1_ && fixed2_) {
            seq_reader<T> r(theta);
            return F::apply(r.read(x1_), x2_);
          } else if (!fixed1_ && !fixed2_) {
            seq_reader<T> r(theta);
            return F::apply(r.read(x1_), r.read(x2_));
          }
          throw std::logic_error("binder_binary illegal state");
        }
      };


      /**
       * Tests whether the binary function specified through the static
       * function F::apply with arguments of the specified types
       * instantiated with double scalars returns right values and
       * first, second, and third derivatives.  It tests all possible
       * combinations of double, var, fvar<double>,
       * fvar<fvar<double>>, fvar<var>, and fvar<fvar<var>>
       * instantiations.
       */
      // test a single sequence of arguments and result
      template <typename F, typename T1, typename T2>
      void test_ad(const T1& x1, const T2& x2, double fx) {
        // create binder then test all autodiff/double combos
        binder_binary<F, T1, T2> f(x1, x2);

        // test (double, double) instantiation
        f.fixed1_ = true;
        f.fixed2_ = true;
        seq_writer<double> a;
        Eigen::VectorXd aaa(0);
        test_functor(f, a.vector(), fx);

        // test (double, autodiff) instantiation
        f.fixed1_ = true;
        f.fixed2_ = false;
        seq_writer<double> b;
        b.write(x2);
        test_functor(f, b.vector(), fx);

        // test (autodiff, double) instantiation
        f.fixed1_ = false;
        f.fixed2_ = true;
        seq_writer<double> c;
        c.write(x1);
        test_functor(f, c.vector(), fx);

        // test (autodiff, autodiff) instantiation
        f.fixed1_ = false;
        f.fixed2_ = false;
        seq_writer<double> d;
        d.write(x1);
        d.write(x2);
        test_functor(f, d.vector(), fx);
      }

    }
  }
}

#endif
