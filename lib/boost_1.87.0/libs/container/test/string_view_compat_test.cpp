//////////////////////////////////////////////////////////////////////////////
//
// (C) Copyright Ion Gaztanaga 2007-2017. Distributed under the Boost
// Software License, Version 1.0. (See accompanying file
// LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/container for documentation.
//
//////////////////////////////////////////////////////////////////////////////

#include <boost/container/string.hpp>
#include <boost/utility/string_view.hpp>

#include <boost/core/lightweight_test.hpp>

#if defined(__cpp_lib_string_view) && (__cpp_lib_string_view >= 201606L)
#define BOOST_CONTAINER_TEST_HAS_STD_STRING_VIEW
#include <string_view>
#endif

template<class StringViewType>
void conversion_test()
{
   #ifndef BOOST_CONTAINER_TEMPLATED_CONVERSION_OPERATOR_BROKEN
   {
      const boost::container::string s = "some text";
      StringViewType sv(s);
      BOOST_TEST(s.data() == sv.data() && s.size() == sv.size());
      StringViewType sv2;
      sv2 = s;
      BOOST_TEST(s.data() == sv2.data() && s.size() == sv2.size());
      const StringViewType csv(s);
      BOOST_TEST(s.data() == sv.data() && s.size() == csv.size());
   }
   #endif
}

template<class StringViewType>
void to_view_test()
{
   const boost::container::string s = "some text";
   StringViewType sv(s.to_view<StringViewType>());
   BOOST_TEST(s.data() == sv.data() && s.size() == sv.size());
   StringViewType sv2;
   sv2 = s.to_view<StringViewType>();
   BOOST_TEST(s.data() == sv2.data() && s.size() == sv2.size());
   const StringViewType csv(s.to_view<StringViewType>());
   BOOST_TEST(s.data() == csv.data() && s.size() == csv.size());
}

template<class StringViewType>
void equal_test()
{
   const StringViewType sv = "same text";
   const StringViewType svd = "different text";
   const boost::container::string s = "same text";
   BOOST_TEST(sv == s);
   BOOST_TEST(s == sv);
   BOOST_TEST(!(svd == s));
   BOOST_TEST(!(s == svd));
}

template<class StringViewType>
void unequal_test()
{
   const StringViewType sv = "same text";
   const StringViewType svd = "different text";
   const boost::container::string s = "same text";
   BOOST_TEST(!(sv != s));
   BOOST_TEST(!(s != sv));
   BOOST_TEST(svd != s);
   BOOST_TEST(s != svd);
}

template<class StringViewType>
void less_test()
{
   StringViewType sv  = "0123456";
   boost::container::string s = "0123459";
   BOOST_TEST(sv < s);
   BOOST_TEST(!(s < sv));

   sv = "0123459";
   s  = "0123456";
   BOOST_TEST(!(sv < s));
   BOOST_TEST(s < sv);

   sv = "0123456";
   BOOST_TEST(!(sv < s));
   BOOST_TEST(!(s < sv));
}

template<class StringViewType>
void greater_test()
{
   StringViewType sv  = "0123459";
   boost::container::string s = "0123456";
   BOOST_TEST(sv > s);
   BOOST_TEST(!(s > sv));

   sv = "0123456";
   s  = "0123459";
   BOOST_TEST(!(sv > s));
   BOOST_TEST(s > sv);

   sv = "0123459";
   BOOST_TEST(!(sv > s));
   BOOST_TEST(!(s > sv));
}

template<class StringViewType>
void less_equal_test()
{
   StringViewType sv  = "0123456";
   boost::container::string s = "0123459";
   BOOST_TEST(sv <= s);
   BOOST_TEST(!(s <= sv));

   sv = "0123459";
   s  = "0123456";
   BOOST_TEST(!(sv <= s));
   BOOST_TEST(s <= sv);

   sv = "0123456";
   BOOST_TEST(sv <= s);
   BOOST_TEST(s <= sv);
}

template<class StringViewType>
void greater_equal_test()
{
   StringViewType sv  = "0123459";
   boost::container::string s = "0123456";
   BOOST_TEST(sv >= s);
   BOOST_TEST(!(s >= sv));

   sv = "0123456";
   s  = "0123459";
   BOOST_TEST(!(sv >= s));
   BOOST_TEST(s >= sv);

   sv = "0123459";
   BOOST_TEST(sv >= s);
   BOOST_TEST(s >= sv);
}

template<class StringViewType>
void constructor_test()
{
   StringViewType sv  = "0123459";
   boost::container::string s(sv);
   BOOST_TEST(sv == s);
   boost::container::string s2(sv, s.get_allocator());
   BOOST_TEST(sv == s);
}

template<class StringViewType>
void assignment_test()
{
   StringViewType sv  = "0123459";
   boost::container::string s;
   s = sv;
   BOOST_TEST(sv == s);
}

template<class StringViewType>
void assign_test()
{
   StringViewType sv  = "0123459";
   boost::container::string s;
   s.assign(sv);
   BOOST_TEST(sv == s);
}

template<class StringViewType>
void plus_equal_test()
{
   StringViewType sv  = "23459";
   boost::container::string s("01");
   s += sv;
   BOOST_TEST(s == "0123459");
}

template<class StringViewType>
void append_test()
{
   StringViewType sv  = "23459";
   boost::container::string s("01");
   s.append(sv);
   BOOST_TEST(s == "0123459");
}

template<class StringViewType>
void insert_test()
{
   StringViewType sv  = "34";
   boost::container::string s("01259");
   s.insert(3u, sv);
   BOOST_TEST(s == "0123459");
}

template<class StringViewType>
void replace_test()
{
   StringViewType sv  = "5678";
   boost::container::string s("01259");
   s.replace(2u, 2u, sv);
   BOOST_TEST(s == "0156789");
   s.replace(s.begin()+3, s.begin()+6, sv);
   BOOST_TEST(s == "01556789");
   s.replace(5u, 3u, sv, 2u, 2u);
   BOOST_TEST(s == "0155678");
}

template<class StringViewType>
void find_test()
{
   const StringViewType sv  = "25";
   boost::container::string s("0125925123");
   BOOST_TEST(s.find(sv,4) == 5);
}

template<class StringViewType>
void rfind_test()
{
   const StringViewType sv  = "25";
   boost::container::string s("0125925123");
   BOOST_TEST(s.rfind(sv,4) == 2);
}

template<class StringViewType>
void find_first_of_test()
{
   const StringViewType sv  = "52";
   boost::container::string s("0125925123");
   BOOST_TEST(s.find_first_of(sv,4) == 5);
}

template<class StringViewType>
void find_last_of_test()
{
   const StringViewType sv  = "52";
   boost::container::string s("520125925123");
   BOOST_TEST(s.find_last_of(sv,6) == 5);
}

template<class StringViewType>
void find_first_not_of_test()
{
   const StringViewType sv  = "52";
   boost::container::string s("0125925123");
   BOOST_TEST(s.find_first_not_of(sv,2) == 4);
}

template<class StringViewType>
void find_last_not_of_test()
{
   const StringViewType sv  = "52";
   boost::container::string s("0125925123");
   BOOST_TEST(s.find_last_not_of(sv,6) == 4);
}

template<class StringViewType>
void compare_test()
{
   const StringViewType sv  = "52";
   boost::container::string s("0125925123");
   BOOST_TEST(s.compare(sv) < 0);
   BOOST_TEST(s.compare(StringViewType("0125925123")) == 0);
   BOOST_TEST(s.compare(2u, s.size() - 2u, StringViewType("25925123")) == 0);
   StringViewType sv2("5212592512389");
   BOOST_TEST(s.compare(2u, s.size() - 2u, sv2, 3, sv2.size()-5u) == 0);
}

#if BOOST_CXX_VERSION >= 201103L
#include <boost/container_hash/hash.hpp>

template<class StringViewType>
void hash_test()
{

   const boost::container::string s1 = "0125925123";
   const boost::container::string s2 = "25925123";
   typedef boost::hash<boost::container::string> shash_t;
   typedef boost::hash<StringViewType> svhash_t;

   BOOST_TEST(shash_t()(s1) == svhash_t()(StringViewType(s1)));
   BOOST_TEST(shash_t()(s2) == svhash_t()(StringViewType(s2)));
}

#else

template<class StringViewType>
void hash_test()
{}

#endif

template<class StringViewType>
void test_all()
{
   conversion_test<StringViewType>();
   to_view_test<StringViewType>();
   equal_test<StringViewType>();
   unequal_test<StringViewType>();
   less_test<StringViewType>();
   greater_test<StringViewType>();
   less_equal_test<StringViewType>();
   greater_equal_test<StringViewType>();
   constructor_test<StringViewType>();
   assignment_test<StringViewType>();
   assign_test<StringViewType>();
   plus_equal_test<StringViewType>();
   append_test<StringViewType>();
   insert_test<StringViewType>();
   replace_test<StringViewType>();
   find_test<StringViewType>();
   rfind_test<StringViewType>();
   find_first_of_test<StringViewType>();
   find_last_of_test<StringViewType>();
   find_first_not_of_test<StringViewType>();
   find_last_not_of_test<StringViewType>();
   compare_test<StringViewType>();
   hash_test<StringViewType>();
}

int main()
{
   test_all<boost::string_view>();
   #ifdef BOOST_CONTAINER_TEST_HAS_STD_STRING_VIEW
   test_all<std::string_view>();
   #endif

   return boost::report_errors();
}

