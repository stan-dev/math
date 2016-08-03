#ifndef STAN_MATH_PRIM_MAT_META_CONTAINER_VIEW_HPP
#define STAN_MATH_PRIM_MAT_META_CONTAINER_VIEW_HPP

#include <stan/math/prim/scal/meta/container_view.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>

namespace stan {
  namespace math {

    /**
     * Template specialization for Eigen::Map view of
     * array with scalar type T2 with size inferred from
     * input Eigen::Matrix    
     *
     * @tparam T1 scalar type of input matrix
     * @tparam T2 scalar type of view.
     * @tparam R rows of input matrix and view
     * @tparam C columns of input matrix and view
     */
    template <typename T1, typename T2, int R, int C>
    class container_view<Eigen::Matrix<T1, R, C>, Eigen::Matrix<T2, R, C> > {
      public:
        /**
         * Initialize Map dimensions with input matrix
         * dimensions
         *
         * @param x input matrix
         * @param y underlying array 
         */
        container_view(const Eigen::Matrix<T1, R, C>& x, T2* y)
         : y_(y, x.rows(), x.cols()) { }

        /**
         * operator[](int i) returns Eigen::Map y
         *
         * @param i index
         */
        Eigen::Map<Eigen::Matrix<T2, R, C> >& operator[](int i) {
          return y_;
        }
      private:
        Eigen::Map<Eigen::Matrix<T2, R, C> > y_;
    };

    /**
     * Template specialization for scalar view of
     * array y with scalar type T2
     *
     * @tparam T1 scalar type of input matrix
     * @tparam T2 scalar type returned by view.
     * @tparam R rows of input matrix and view
     * @tparam C columns of input matrix and view
     */

    template <typename T1, typename T2, int R, int C>
    class container_view<Eigen::Matrix<T1, R, C>, T2> {
      public:
        /** 
         * Constructor
         *
         * @param x input matrix
         * @param y underlying array 
         */
        container_view(const Eigen::Matrix<T1, R, C>& x, T2* y)
         : y_(y) { }

        /**
         * operator[](int i) returns reference to scalar 
         * of type T2 at appropriate index i in array y
        */
        T2& operator[](int i) {
          return y_[i];
        }
      private:
        T2* y_;
    };

    /**
     * Template specialization for matrix view of
     * array y with scalar type T2 with shape
     * equal to x
     *
     * @tparam T1 scalar type of input vector of matrices
     * @tparam T2 scalar type of matrix view
     * @tparam R rows of input matrix and view
     * @tparam C columns of input matrix and view
     */
    template <typename T1, typename T2, int R, int C>
    class container_view<std::vector<Eigen::Matrix<T1, R, C> >,
                          Eigen::Matrix<T2, R, C> > {
      public:
        /** 
         * Constructor assumes all matrix elements in 
         * std::vector are of same dimension
         *
         * Initializes y_view as 1x1 matrix because no
         * nullary constructor for Eigen::Map
         *
         * @param x input matrix
         * @param y underlying array 
         */
        container_view(const std::vector<Eigen::Matrix<T1, R, C> >& x, T2* y)
         : y_view(y, 1, 1), y_(y) {
           if (x.size() > 0) {
             rows = x[0].rows();
             cols = x[0].cols();
           } else {
             rows = 0;
             cols = 0;
           }
         }

        /**
         * operator[](int i) returns matrix view  
         * of scalartype T2 at appropriate index i in array y
        */
        Eigen::Map<Eigen::Matrix<T2, R, C> >& operator[](int i) {
          int offset = i * rows * cols;
          new (&y_view) Eigen::Map<Eigen::Matrix<T2, R, C> >
            (y_ + offset, rows, cols);
          return y_view;
        }
      private:
        Eigen::Map<Eigen::Matrix<T2, R, C> > y_view;
        T2* y_;
        int rows;
        int cols;
    };
  }
}
#endif
