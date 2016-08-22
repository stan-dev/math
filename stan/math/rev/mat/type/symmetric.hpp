#ifndef STAN_MATH_REV_MAT_TYPE_SYMMETRIC_HPP
#define STAN_MATH_REV_MAT_TYPE_SYMMETRIC_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/scal/fun/value_of.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>

namespace stan {
  namespace math {

    /**
     * Symmetric matrix class.
     *
     * Stores vars for a symmetric matrix.
     * Only stores lower-triangular values.
     *
     * @param matrix A
     * @param matrix L, cholesky factor of A
     */
    class symmetric {
      private:
        int rows_;
        int size_;

        inline int vech_ind(int i, int j, int N) {
          int ind = j * N + i - (j * (j + 1)) / 2;
          return ind;
        }
      public:
        Eigen::Matrix<var, -1, 1> vars_;

        symmetric(const Eigen::Matrix<var, -1, -1>& vars,
                  const int rows)
          : rows_(rows), size_(rows_ * (rows_ + 1) / 2) {
            check_size_match("symmetric",
                             "Elements of ", "data", vars.size(),
                             "value ", "rows * (rows + 1) / 2", size_);
            vars_.resize(size_);
            for (int i = 0; i < size_; ++i) 
              vars_.coeffRef(i) = vars.coeffRef(i);
          }

        symmetric(const int rows)
          : rows_(rows), size_(rows_ * (rows_ + 1) / 2) {
            vars_.resize(size_);
          } 

        var& operator()(int i, int j) {
          check_less_or_equal("Indexer", "rows", j, i);
          int ind = vech_ind(i, j, rows_);
          check_less_or_equal("Indexer", "element", ind, size_);
          return vars_.coeffRef(ind);
        }

        int size() {
          return size_;
        }
    };
  }
}
#endif
