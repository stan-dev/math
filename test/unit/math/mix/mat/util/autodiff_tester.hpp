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

      template <typename T1, typename T2>
      void expect_near(const T1& x1, const T2& x2,
                       double tol = 1e-11) {
        if (is_nan(x1) || is_nan(x2))
          EXPECT_TRUE(is_nan(x1) && is_nan(x2));
        else if (is_inf(x1) || is_inf(x2))
          EXPECT_EQ(x1, x2);
        else
          EXPECT_NEAR(x1, x2, tol);
      }


      template <typename T, int R, int C>
      void expect_near(const Eigen::Matrix<T, R, C>& x1,
                       const Eigen::Matrix<T, R, C>& x2,
                       double tol = 1e-9) {
        EXPECT_EQ(x1.rows(), x2.rows());
        EXPECT_EQ(x1.cols(), x2.cols());
        for (int i = 0; i < x1.size(); ++i)
          expect_near(x1(i), x2(i), tol);
      }


      template <typename F>
      void test_value(const F& f, const Eigen::VectorXd& x, double fx) {
        if (is_nan(fx))
          EXPECT_TRUE(is_nan(f(x)));
        else
          EXPECT_FLOAT_EQ(fx, f(x));
      }


      // var
      template <typename F>
      void test_gradient(const F& f, const Eigen::VectorXd& x, double fx) {
        Eigen::VectorXd grad_ad;
        double fx_ad;
        gradient(f, x, fx_ad, grad_ad);
        expect_near(fx, fx_ad);
        Eigen::VectorXd grad_fd;
        double fx_fd;
        finite_diff_gradient(f, x, fx_fd, grad_fd);
        expect_near(grad_fd, grad_ad);
      }


      // var
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


      // fvar<fvar<double>>
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


      // fvar<var>
      template <typename F>
      void test_hessian(const F& f, const Eigen::VectorXd& x, double fx) {
        double fx_ad;
        Eigen::VectorXd grad_ad;
        Eigen::MatrixXd H_ad;
        hessian(f, x, fx_ad, grad_ad, H_ad);
        expect_near(fx, fx_ad);
        double fx_fd;
        Eigen::VectorXd grad_fd;
        Eigen::MatrixXd H_fd;
        finite_diff_hessian(f, x, fx_fd, grad_fd, H_fd);
        expect_near(grad_fd, grad_ad);
        expect_near(H_fd, H_ad);
      }


      // fvar<fvar<var>>
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
        test_grad_hessian(f, x, fx);
      }


      // adapts two-argument function to vector function
      // with x1 and x2 fixed or as arguments in the vector
      template <typename F, typename T1, typename T2>
      struct binder_binary {
        // stored fixed values
        T1 x1_;
        T2 x2_;

        // flags indicating whether to use fixed values
        bool fixed1_;
        bool fixed2_;

        binder_binary(const T1& x1, const T2& x2)
          : x1_(x1), x2_(x2), fixed1_(false), fixed2_(false) { }

        template <typename T>
        T operator()(const Eigen::Matrix<T, -1, 1>& theta) const {
          if (fixed1_ && !fixed2_) {
            seq_reader<T> r(theta);
            return F::apply(x1_, r.read(x2_));
          } else if (!fixed1_ && fixed2_) {
            seq_reader<T> r(theta);
            return F::apply(r.read(x1_), x2_);
          } else if (!fixed1_ && !fixed2_) {
            seq_reader<T> r(theta);
            return F::apply(r.read(x1_), r.read(x2_));
          }
          throw std::logic_error("fixed1_ = true and fixed2_ = true illegal");
          return T();
        }
      };


      // test a single sequence of arguments and result
      template <typename F, typename T1, typename T2>
      void test_ad(const T1& x1, const T2& x2, double fx) {
        // create binder then test all autodiff/double combos
        binder_binary<F, T1, T2> f(x1, x2);

        // double, autodiff
        f.fixed1_ = true;
        f.fixed2_ = false;
        seq_writer<double> a;
        a.write(x2);
        test_functor(f, a.vector(), fx);

        // autodiff, double
        f.fixed1_ = false;
        f.fixed2_ = true;
        seq_writer<double> b;
        b.write(x1);
        test_functor(f, b.vector(), fx);

        // autodiff, autodiff
        f.fixed1_ = false;
        f.fixed2_ = false;
        seq_writer<double> c;
        c.write(x1);
        c.write(x2);
        test_functor(f, c.vector(), fx);
      }

    }
  }
}

#endif
