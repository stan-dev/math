#pragma once


namespace boost {

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/release/libs/assert/doc/html/assert.html#source_location
struct source_location {};

namespace container {
namespace pmr {

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/release/doc/html/boost/container/pmr/polymorphic_allocator.html
template <class T>
struct polymorphic_allocator {};

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/1_85_0/doc/html/boost/container/pmr/memory_resource.html
struct memory_resource {};

} // namespace pmr
} // namespace container

namespace system {

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/release/libs/system/doc/html/system.html#ref_error_code
struct error_code {};

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/release/libs/system/doc/html/system.html#ref_error_category
struct error_category {};

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/release/libs/system/doc/html/system.html#ref_error_condition
struct error_condition {};

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/release/libs/system/doc/html/system.html#ref_system_error
struct system_error {};

/// !EXTERNAL!
///
/// @see https://www.boost.org/doc/libs/release/libs/system/doc/html/system.html#ref_boostsystemresult_hpp
struct result {};

} // namespace system
} // namespace boost

namespace std {

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/iterator/reverse_iterator
template <class T>
struct reverse_iterator {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/types/size_t
struct size_t {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/types/integer
struct uint64_t {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/types/integer
struct int64_t {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/types/nullptr_t
struct nullptr_t {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/types/ptrdiff_t
struct ptrdiff_t {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/utility/initializer_list
template <class T>
struct initializer_list {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/error/error_code
struct error_code {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/utility/pair
struct pair {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/types/byte
struct byte {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/io/basic_ostream
struct ostream {};

/// !EXTERNAL!
///
/// @see https://en.cppreference.com/w/cpp/io/basic_istream
struct istream {};

} // namespace std
