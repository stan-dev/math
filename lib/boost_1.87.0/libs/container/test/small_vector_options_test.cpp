//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2004-2013. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/container for documentation.
//
//////////////////////////////////////////////////////////////////////////////
#include <boost/container/small_vector.hpp>
#include <boost/container/allocator.hpp>
#include <boost/core/lightweight_test.hpp>
#include <boost/assert.hpp>
using namespace boost::container;

const std::size_t Capacity = 10u;

void test_alignment()
{
   {  //extended alignment
      const std::size_t extended_alignment = sizeof(int)*4u;
      BOOST_CONTAINER_STATIC_ASSERT(extended_alignment > dtl::alignment_of<int>::value);
      #if !defined(BOOST_NO_CXX11_TEMPLATE_ALIASES)
      using options_t = small_vector_options_t< inplace_alignment<extended_alignment> >;
      #else
      typedef small_vector_options
         < inplace_alignment<extended_alignment> >::type options_t;
      #endif

      small_vector<int, Capacity, void, options_t> v;
      v.resize(v.capacity());
      BOOST_ASSERT((reinterpret_cast<std::size_t>(&v[0]) % extended_alignment) == 0);
   }
   {  //default alignment
      #if !defined(BOOST_NO_CXX11_TEMPLATE_ALIASES)
      using options_t = small_vector_options_t< inplace_alignment<0> >;
      #else
      typedef small_vector_options< inplace_alignment<0> >::type options_t;
      #endif

      small_vector<int, Capacity, void, options_t> v;
      v.resize(v.capacity());
      BOOST_ASSERT((reinterpret_cast<std::size_t>(&v[0]) % dtl::alignment_of<int>::value) == 0);
   }
}

void test_growth_factor_50()
{
   #if !defined(BOOST_NO_CXX11_TEMPLATE_ALIASES)
   using options_t = small_vector_options_t< growth_factor<growth_factor_50> >;
   #else
   typedef small_vector_options
      < growth_factor<growth_factor_50> >::type options_t;
   #endif

   small_vector<int, Capacity, new_allocator<int>, options_t> v;

   v.resize(5);
   v.resize(v.capacity());
   std::size_t old_capacity = v.capacity();
   v.push_back(0);
   std::size_t new_capacity = v.capacity();
   BOOST_TEST(new_capacity == old_capacity + old_capacity/2);
}

void test_growth_factor_60()
{
   #if !defined(BOOST_NO_CXX11_TEMPLATE_ALIASES)
   using options_t = small_vector_options_t< growth_factor<growth_factor_60> >;
   #else
   typedef small_vector_options
      < growth_factor<growth_factor_60> >::type options_t;
   #endif

   small_vector<int, Capacity, new_allocator<int>, options_t> v;

   v.resize(5);
   v.resize(v.capacity());
   std::size_t old_capacity = v.capacity();
   v.push_back(0);
   std::size_t new_capacity = v.capacity();
   BOOST_TEST(new_capacity == old_capacity + 3*old_capacity/5);
}

void test_growth_factor_100()
{
   #if !defined(BOOST_NO_CXX11_TEMPLATE_ALIASES)
   using options_t = small_vector_options_t< growth_factor<growth_factor_100> >;
   #else
   typedef small_vector_options
      < growth_factor<growth_factor_100> >::type options_t;
   #endif

   small_vector<int, Capacity, new_allocator<int>, options_t> v;

   v.resize(5);
   v.resize(v.capacity());
   std::size_t old_capacity = v.capacity();
   v.push_back(0);
   std::size_t new_capacity = v.capacity();
   BOOST_TEST(new_capacity == 2*old_capacity);
}

template<class Unsigned, class VectorType>
void test_stored_size_type_impl()
{
   #ifndef BOOST_NO_EXCEPTIONS
   VectorType v;
   typedef typename VectorType::size_type    size_type;
   typedef typename VectorType::value_type   value_type;
   size_type const max = Unsigned(-1);
   v.resize(5);
   v.resize(max);
   BOOST_TEST_THROWS(v.resize(max+1),                    std::exception);
   BOOST_TEST_THROWS(v.push_back(value_type(1)),         std::exception);
   BOOST_TEST_THROWS(v.insert(v.begin(), value_type(1)), std::exception);
   BOOST_TEST_THROWS(v.emplace(v.begin(), value_type(1)),std::exception);
   BOOST_TEST_THROWS(v.reserve(max+1),                   std::exception);
   BOOST_TEST_THROWS(VectorType v2(max+1),               std::exception);
   #endif
}

template<class Unsigned>
void test_stored_size_type()
{
   #if !defined(BOOST_NO_CXX11_TEMPLATE_ALIASES)
   using options_t = small_vector_options_t< stored_size<Unsigned> >;
   #else
   typedef typename small_vector_options
      < stored_size<Unsigned> >::type options_t;
   #endif

   typedef small_vector<unsigned char, Unsigned(-1)> normal_small_vector_t;

   //Test first with a typical allocator
   {
      typedef small_vector<unsigned char, Unsigned(-1), new_allocator<unsigned char>, options_t> small_vector_t;
      test_stored_size_type_impl<Unsigned, small_vector_t>();
      BOOST_CONTAINER_STATIC_ASSERT(sizeof(normal_small_vector_t) > sizeof(small_vector_t));
   }
   //Test with a V2 allocator
   {
      typedef small_vector<unsigned char, Unsigned(-1), allocator<unsigned char>, options_t> small_vector_t;
      test_stored_size_type_impl<Unsigned, small_vector_t>();
      BOOST_CONTAINER_STATIC_ASSERT(sizeof(normal_small_vector_t) > sizeof(small_vector_t));
   }
}

int main()
{
   test_alignment();
   test_growth_factor_50();
   test_growth_factor_60();
   test_growth_factor_100();
   test_stored_size_type<unsigned char>();
   test_stored_size_type<unsigned short>();
   return ::boost::report_errors();
}
