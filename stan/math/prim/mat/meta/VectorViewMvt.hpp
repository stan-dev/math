#ifndef STAN_MATH_PRIM_MAT_META_VECTORVIEWMVT_HPP
#define STAN_MATH_PRIM_MAT_META_VECTORVIEWMVT_HPP

#include <stan/math/prim/mat/meta/is_vector_like.hpp>
#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/scal/meta/scalar_type_pre.hpp>
#include <stdexcept>
#include <vector>

namespace stan {


  /**
   * VectorViewMvt is a template expression that wraps either an
   * Eigen::Matrix or a std::vector<Eigen::Matrix> and allows the
   * template expression to be used as an array using
   * <code>operator[]</code>.
   *
   * @tparam T Type of scalar of the matrix being wrapped.  
   * @tparam is_array True if underlying type T can be indexed with
   *   operator[].  
   * @tparam throw_if_accessed True if the behavior is to
   *   throw an exception whenever <code>operator[]</code> is called.
   */
  template <typename T, bool is_array
            = stan::is_vector_like
            <typename stan::math::value_type<T>::type>::value,
            bool throw_if_accessed = false>
  class VectorViewMvt {
  public:
    typedef typename scalar_type_pre<T>::type matrix_t;

    /**
     * Constructor.
     */
    explicit VectorViewMvt(matrix_t& m) : x_(&m) { }

    /**
     * Constructor.
     */
    explicit VectorViewMvt(std::vector<matrix_t>& vm) : x_(&vm[0]) { }

    /**
     * Allows the structure to be accessed like an array. If is_array
     * is false, this will return the matrix it was constructed
     * with. If is_array is true, This does not check bounds and will
     * likely segfault if the index is out of range.
     *
     * @param i index. Only used if access is true.
     * @return Reference to a matrix.
     * @throw std::out_of_range if the template parameter, 
     *   throw_if_accessed, is true.
     */
    matrix_t& operator[](int i) {
      if (throw_if_accessed)
        throw std::out_of_range("VectorViewMvt: this cannot be accessed");
      if (is_array)
        return x_[i];
      else
        return x_[0];
    }

  private:
    matrix_t* x_;
  };

  /**
   * VectorViewMvt with const correctness.
   *
   */
  template <typename T, bool is_array, bool throw_if_accessed>
  class VectorViewMvt<const T, is_array, throw_if_accessed> {
  public:
    typedef typename scalar_type_pre<T>::type matrix_t;

    explicit VectorViewMvt(const matrix_t& m) : x_(&m) { }

    explicit VectorViewMvt(const std::vector<matrix_t>& vm) : x_(&vm[0]) { }

    /**
     * Allows the structure to be accessed like an array. If is_array
     * is false, this will return the matrix it was constructed
     * with. If is_array is true, This does not check bounds and will
     * likely segfault if the index is out of range.
     *
     * @param i index. Only used if access is true.
     * @return Reference to a matrix.
     * @throw std::out_of_range if the template parameter, 
     *   throw_if_accessed, is true.
     */
    const matrix_t& operator[](int i) const {
      if (throw_if_accessed)
        throw std::out_of_range("VectorViewMvt: this cannot be accessed");
      if (is_array)
        return x_[i];
      else
        return x_[0];
    }

  private:
    const matrix_t* x_;
  };

}
#endif

