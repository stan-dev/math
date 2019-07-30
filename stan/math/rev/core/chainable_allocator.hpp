#ifndef STAN_MATH_REV_CORE_CHAINABLE_ALLOCATOR_HPP
#define STAN_MATH_REV_CORE_CHAINABLE_ALLOCATOR_HPP

#include <stan/math/rev/core/chainablestack.hpp>
#include <memory>

namespace stan {
namespace math {

/** \class chainable_allocator
* \ingroup Core_Module
*
* \brief STL compatible allocator to use with types requiring a non standrad alignment.
*
* The memory is aligned as for dynamically aligned matrix/array types such as MatrixXd.
* By default, it will thus provide at least 16 bytes alignment and more in following cases:
*  - 32 bytes alignment if AVX is enabled.
*  - 64 bytes alignment if AVX512 is enabled.
*
* This can be controlled using the \c EIGEN_MAX_ALIGN_BYTES macro as documented
* \link TopicPreprocessorDirectivesPerformance there \endlink.
*
* Example:
* \code
* // Matrix4f requires 16 bytes alignment:
* std::map< int, Matrix4f, std::less<int>,
*           chainable_allocator<std::pair<const int, Matrix4f> > > my_map_mat4;
* // Vector3f does not require 16 bytes alignment, no need to use Eigen's allocator:
* std::map< int, Vector3f > my_map_vec3;
* \endcode
*
* \sa \blank \ref TopicStlContainers.
*/
template<class T>
class chainable_allocator : public std::allocator<T> {
public:
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;
  typedef T*              pointer;
  typedef const T*        const_pointer;
  typedef T&              reference;
  typedef const T&        const_reference;
  typedef T               value_type;

  template<class U>
  struct rebind
  {
    typedef chainable_allocator<U> other;
  };

  chainable_allocator() : std::allocator<T>() {}

  chainable_allocator(const chainable_allocator& other) : std::allocator<T>(other) {}

  template<class U>
  chainable_allocator(const chainable_allocator<U>& other) : std::allocator<T>(other) {}

  ~chainable_allocator() {}

  // not copyable or movable
  chainable_allocator& operator=(const chainable_allocator&) = delete;
  chainable_allocator& operator=(chainable_allocator&&) = delete;


  pointer allocate(size_type num, const void* /*hint*/ = 0)
  {
    return static_cast<pointer>(ChainableStack::instance_->memalloc_.alloc(num * sizeof(T)));
  }

  void deallocate(pointer p, size_type /*num*/)
  {
    //Nope!
  }
};

}  // namespace math
}  // namespace stan
#endif
