#ifndef STAN_MATH_REV_CORE_CHAINABLESTACK_HPP
#define STAN_MATH_REV_CORE_CHAINABLESTACK_HPP

#include <stan/math/rev/core/autodiffstackstorage.hpp>
namespace stan {
namespace math {
class chainable_alloc;
class vari_base;
template <typename T>
class var_value;
using ChainableStack = AutodiffStackSingleton<vari_base, chainable_alloc>;

}  // namespace math
}  // namespace stan

namespace Eigen {
namespace internal {
/** \internal Allocates \a size objects of type T. The returned pointer is
 * guaranteed to have 16 bytes alignment. On allocation error, the returned
 * pointer is undefined, but a std::bad_alloc is thrown. The default constructor
 * of T is called.
 */
EIGEN_DEVICE_FUNC inline stan::math::var_value<double> *aligned_new(
    std::size_t size) {
  using T = stan::math::var_value<double>;
  check_size_for_overflow<T>(size);
  T *result = reinterpret_cast<T *>(aligned_malloc(sizeof(T) * size));
  return result;
}

template <>
EIGEN_DEVICE_FUNC inline stan::math::var_value<double> *
conditional_aligned_new<stan::math::var_value<double>, true>(std::size_t size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = true;
  check_size_for_overflow<T>(size);
  T *result = reinterpret_cast<T *>(
      conditional_aligned_malloc<Align>(sizeof(T) * size));
  EIGEN_TRY { return construct_elements_of_array(result, size); }
  EIGEN_CATCH(...) {
    conditional_aligned_free<Align>(result);
    EIGEN_THROW;
  }
  return result;
}

template <>
EIGEN_DEVICE_FUNC inline stan::math::var_value<double>
    *conditional_aligned_new<stan::math::var_value<double>, false>(
        std::size_t size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = false;
  check_size_for_overflow<T>(size);
  T *result = reinterpret_cast<T *>(
      conditional_aligned_malloc<Align>(sizeof(T) * size));
  EIGEN_TRY { return construct_elements_of_array(result, size); }
  EIGEN_CATCH(...) {
    conditional_aligned_free<Align>(result);
    EIGEN_THROW;
  }
  return result;
}

/** \internal Deletes objects constructed with aligned_new
 * The \a size parameters tells on how many objects to call the destructor of T.
 */
EIGEN_DEVICE_FUNC inline void aligned_delete(stan::math::var_value<double> *ptr,
                                             std::size_t size) {
  using T = stan::math::var_value<double>;
  destruct_elements_of_array<T>(ptr, size);
  aligned_free(ptr);
}

/** \internal Deletes objects constructed with conditional_aligned_new
 * The \a size parameters tells on how many objects to call the destructor of T.
 */
template <>
EIGEN_DEVICE_FUNC inline void
conditional_aligned_delete<stan::math::var_value<double>, true>(
    stan::math::var_value<double> *ptr, std::size_t size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = true;
  destruct_elements_of_array<T>(ptr, size);
  conditional_aligned_free<Align>(ptr);
}

template <>
EIGEN_DEVICE_FUNC inline void
conditional_aligned_delete<stan::math::var_value<double>, false>(
    stan::math::var_value<double> *ptr, std::size_t size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = false;
  destruct_elements_of_array<T>(ptr, size);
  conditional_aligned_free<Align>(ptr);
}

template <>
EIGEN_DEVICE_FUNC inline stan::math::var_value<double>
    *conditional_aligned_realloc_new<stan::math::var_value<double>, true>(
        stan::math::var_value<double> *pts, std::size_t new_size,
        std::size_t old_size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = true;
  check_size_for_overflow<T>(new_size);
  check_size_for_overflow<T>(old_size);
  if (new_size < old_size)
    destruct_elements_of_array(pts + new_size, old_size - new_size);
  T *result = reinterpret_cast<T *>(conditional_aligned_realloc<Align>(
      reinterpret_cast<void *>(pts), sizeof(T) * new_size,
      sizeof(T) * old_size));
  if (new_size > old_size) {
    EIGEN_TRY {
      construct_elements_of_array(result + old_size, new_size - old_size);
    }
    EIGEN_CATCH(...) {
      conditional_aligned_free<Align>(result);
      EIGEN_THROW;
    }
  }
  return result;
}

template <>
EIGEN_DEVICE_FUNC inline stan::math::var_value<double>
    *conditional_aligned_realloc_new<stan::math::var_value<double>, false>(
        stan::math::var_value<double> *pts, std::size_t new_size,
        std::size_t old_size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = false;
  check_size_for_overflow<T>(new_size);
  check_size_for_overflow<T>(old_size);
  if (new_size < old_size)
    destruct_elements_of_array(pts + new_size, old_size - new_size);
  T *result = reinterpret_cast<T *>(conditional_aligned_realloc<Align>(
      reinterpret_cast<void *>(pts), sizeof(T) * new_size,
      sizeof(T) * old_size));
  if (new_size > old_size) {
    EIGEN_TRY {
      construct_elements_of_array(result + old_size, new_size - old_size);
    }
    EIGEN_CATCH(...) {
      conditional_aligned_free<Align>(result);
      EIGEN_THROW;
    }
  }
  return result;
}

template <>
EIGEN_DEVICE_FUNC inline stan::math::var_value<double>
    *conditional_aligned_new_auto<stan::math::var_value<double>, true>(
        std::size_t size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = true;
  if (size == 0)
    return 0;  // short-cut. Also fixes Bug 884
  check_size_for_overflow<T>(size);
  T *result = reinterpret_cast<T *>(
      conditional_aligned_malloc<Align>(sizeof(T) * size));
  if (NumTraits<T>::RequireInitialization) {
    EIGEN_TRY { construct_elements_of_array(result, size); }
    EIGEN_CATCH(...) {
      conditional_aligned_free<Align>(result);
      EIGEN_THROW;
    }
  }
  return result;
}

template <>
EIGEN_DEVICE_FUNC inline stan::math::var_value<double>
    *conditional_aligned_new_auto<stan::math::var_value<double>, false>(
        std::size_t size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = false;
  if (size == 0)
    return 0;  // short-cut. Also fixes Bug 884
  check_size_for_overflow<T>(size);
  T *result = reinterpret_cast<T *>(
      conditional_aligned_malloc<Align>(sizeof(T) * size));
  if (NumTraits<T>::RequireInitialization) {
    EIGEN_TRY { construct_elements_of_array(result, size); }
    EIGEN_CATCH(...) {
      conditional_aligned_free<Align>(result);
      EIGEN_THROW;
    }
  }
  return result;
}

template <>
inline stan::math::var_value<double>
    *conditional_aligned_realloc_new_auto<stan::math::var_value<double>, true>(
        stan::math::var_value<double> *pts, std::size_t new_size,
        std::size_t old_size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = true;
  check_size_for_overflow<T>(new_size);
  check_size_for_overflow<T>(old_size);
  if (NumTraits<T>::RequireInitialization && (new_size < old_size))
    destruct_elements_of_array(pts + new_size, old_size - new_size);
  T *result = reinterpret_cast<T *>(conditional_aligned_realloc<Align>(
      reinterpret_cast<void *>(pts), sizeof(T) * new_size,
      sizeof(T) * old_size));
  if (NumTraits<T>::RequireInitialization && (new_size > old_size)) {
    EIGEN_TRY {
      construct_elements_of_array(result + old_size, new_size - old_size);
    }
    EIGEN_CATCH(...) {
      conditional_aligned_free<Align>(result);
      EIGEN_THROW;
    }
  }
  return result;
}

template <>
inline stan::math::var_value<double>
    *conditional_aligned_realloc_new_auto<stan::math::var_value<double>, false>(
        stan::math::var_value<double> *pts, std::size_t new_size,
        std::size_t old_size) {
  using T = stan::math::var_value<double>;
  constexpr bool Align = false;
  check_size_for_overflow<T>(new_size);
  check_size_for_overflow<T>(old_size);
  if (NumTraits<T>::RequireInitialization && (new_size < old_size))
    destruct_elements_of_array(pts + new_size, old_size - new_size);
  T *result = reinterpret_cast<T *>(conditional_aligned_realloc<Align>(
      reinterpret_cast<void *>(pts), sizeof(T) * new_size,
      sizeof(T) * old_size));
  if (NumTraits<T>::RequireInitialization && (new_size > old_size)) {
    EIGEN_TRY {
      construct_elements_of_array(result + old_size, new_size - old_size);
    }
    EIGEN_CATCH(...) {
      conditional_aligned_free<Align>(result);
      EIGEN_THROW;
    }
  }
  return result;
}

template <>
EIGEN_DEVICE_FUNC inline void conditional_aligned_delete_auto<stan::math::var_value<double>, true>(
    stan::math::var_value<double>* ptr, std::size_t size) {
  conditional_aligned_free<true>(ptr);
}

template <>
EIGEN_DEVICE_FUNC inline void conditional_aligned_delete_auto<stan::math::var_value<double>, false>(
    stan::math::var_value<double>* ptr, std::size_t size) {
  conditional_aligned_free<false>(ptr);
}


}  // namespace internal
}  // namespace Eigen
#endif
