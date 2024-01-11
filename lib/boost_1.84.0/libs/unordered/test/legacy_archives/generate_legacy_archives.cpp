/* Copyright 2023 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/unordered for library home page.
 */

/* This program has been used to generate archives of Boost.Unordered
 * containers with Boost 1.83, when serialization support was provided
 * externally to Boost.Unordered in
 * boost/serialization/boost_unordered_(map|set).hpp. Beginning with the
 * release of native Boost.Unordered serialization capabilities in Boost
 * 1.84, these archives are used to test backwards loading compatibility
 * as enabled by BOOST_UNORDERED_ENABLE_SERIALIZATION_COMPATIBILITY_V0.
 */

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/boost_unordered_map.hpp>
#include <boost/serialization/boost_unordered_set.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/version.hpp>
#include <fstream>
#include <random>

template<typename Value=unsigned int>
struct random_generator
{
  Value operator()()
  {
    return static_cast<Value>(dist(gen));
  }

  std::uniform_int_distribution<unsigned int> dist;
  std::mt19937                                gen{231};
};

template<>
struct random_generator<std::string>
{
  std::string operator()()
  {
    return std::to_string(rng());
  }

  random_generator<> rng;
};

template<>
struct random_generator<const std::string>:random_generator<std::string>{};

template<typename T,typename Q>
struct random_generator<std::pair<T,Q>>
{
  std::pair<T,Q> operator()()
  {
    return {rngt(),rngq()};
  }

  random_generator<T> rngt;
  random_generator<Q> rngq;
};

template<typename T>
struct non_const
{
  typedef T type;
};

template<typename T>
struct non_const<const T>
{
  using type=typename non_const<T>::type;
};

template<typename T, typename Q>
struct non_const<std::pair<T, Q>>
{
  using type=std::pair<
    typename non_const<T>::type,
    typename non_const<Q>::type>;
};

template<typename Container,typename Archive>
void generate_legacy_archive(const char* filename,std::size_t n)
{
  using value_type=typename Container::value_type;
  using vector=std::vector<typename non_const<value_type>::type>;

  Container                       c;
  vector                          v;
  random_generator<value_type>    rng;
  std::uniform_int_distribution<> repeat(0,1);
  std::mt19937                    gen{231};

  for(std::size_t i=0;i<n;++i){
    value_type x=rng();
    c.insert(x);
    v.push_back(x);
    if(repeat(gen)){
      c.insert(x);
      v.push_back(x);
    }
  }

  std::ofstream ofs(filename);
  Archive oa(ofs);
  oa<<boost::serialization::make_nvp("container",c);
  oa<<boost::serialization::make_nvp("values",v);
}

int main()
{
  static_assert(BOOST_VERSION<=108300,"to be used with Boost <1.84.");

  generate_legacy_archive<
    boost::unordered_map<int,int>,boost::archive::text_oarchive
  >("map_int_0.txt",0);
  generate_legacy_archive<
    boost::unordered_map<int,int>,boost::archive::text_oarchive
  >("map_int_10.txt",10);
  generate_legacy_archive<
    boost::unordered_map<int,int>,boost::archive::text_oarchive
  >("map_int_100.txt",100);
  generate_legacy_archive<
    boost::unordered_map<std::string,std::string>,boost::archive::text_oarchive
  >("map_string_0.txt",0);
  generate_legacy_archive<
    boost::unordered_map<std::string,std::string>,boost::archive::text_oarchive
  >("map_string_10.txt",10);
  generate_legacy_archive<boost::unordered_map<
    std::string,std::string>,boost::archive::text_oarchive
  >("map_string_100.txt",100);

  generate_legacy_archive<
    boost::unordered_multimap<int,int>,boost::archive::text_oarchive
  >("multimap_int_0.txt",0);
  generate_legacy_archive<
    boost::unordered_multimap<int,int>,boost::archive::text_oarchive
  >("multimap_int_10.txt",10);
  generate_legacy_archive<
    boost::unordered_multimap<int,int>,boost::archive::text_oarchive
  >("multimap_int_100.txt",100);
  generate_legacy_archive<
    boost::unordered_multimap<std::string,std::string>,boost::archive::text_oarchive
  >("multimap_string_0.txt",0);
  generate_legacy_archive<
    boost::unordered_multimap<std::string,std::string>,boost::archive::text_oarchive
  >("multimap_string_10.txt",10);
  generate_legacy_archive<
    boost::unordered_multimap<std::string,std::string>,boost::archive::text_oarchive
  >("multimap_string_100.txt",100);

  generate_legacy_archive<
    boost::unordered_set<int>,boost::archive::text_oarchive
  >("set_int_0.txt",0);
  generate_legacy_archive<
    boost::unordered_set<int>,boost::archive::text_oarchive
  >("set_int_10.txt",10);
  generate_legacy_archive<
    boost::unordered_set<int>,boost::archive::text_oarchive
  >("set_int_100.txt",100);
  generate_legacy_archive<
    boost::unordered_set<std::string>,boost::archive::text_oarchive
  >("set_string_0.txt",0);
  generate_legacy_archive<
    boost::unordered_set<std::string>,boost::archive::text_oarchive
  >("set_string_10.txt",10);
  generate_legacy_archive<
    boost::unordered_set<std::string>,boost::archive::text_oarchive
  >("set_string_100.txt",100);

  generate_legacy_archive<
    boost::unordered_multiset<int>,boost::archive::text_oarchive
  >("multiset_int_0.txt",0);
  generate_legacy_archive<
    boost::unordered_multiset<int>,boost::archive::text_oarchive
  >("multiset_int_10.txt",10);
  generate_legacy_archive<
    boost::unordered_multiset<int>,boost::archive::text_oarchive
  >("multiset_int_100.txt",100);
  generate_legacy_archive<
    boost::unordered_multiset<std::string>,boost::archive::text_oarchive
  >("multiset_string_0.txt",0);
  generate_legacy_archive<
    boost::unordered_multiset<std::string>,boost::archive::text_oarchive
  >("multiset_string_10.txt",10);
  generate_legacy_archive<
    boost::unordered_multiset<std::string>,boost::archive::text_oarchive
  >("multiset_string_100.txt",100);

  generate_legacy_archive<
    boost::unordered_map<int,int>,boost::archive::xml_oarchive
  >("map_int_0.xml",0);
  generate_legacy_archive<
    boost::unordered_map<int,int>,boost::archive::xml_oarchive
  >("map_int_10.xml",10);
  generate_legacy_archive<
    boost::unordered_map<int,int>,boost::archive::xml_oarchive
  >("map_int_100.xml",100);
  generate_legacy_archive<
    boost::unordered_map<std::string,std::string>,boost::archive::xml_oarchive
  >("map_string_0.xml",0);
  generate_legacy_archive<
    boost::unordered_map<std::string,std::string>,boost::archive::xml_oarchive
  >("map_string_10.xml",10);
  generate_legacy_archive<
    boost::unordered_map<std::string,std::string>,boost::archive::xml_oarchive
  >("map_string_100.xml",100);

  generate_legacy_archive<
    boost::unordered_multimap<int,int>,boost::archive::xml_oarchive
  >("multimap_int_0.xml",0);
  generate_legacy_archive<
    boost::unordered_multimap<int,int>,boost::archive::xml_oarchive
  >("multimap_int_10.xml",10);
  generate_legacy_archive<
    boost::unordered_multimap<int,int>,boost::archive::xml_oarchive
  >("multimap_int_100.xml",100);
  generate_legacy_archive<
    boost::unordered_multimap<std::string,std::string>,boost::archive::xml_oarchive
  >("multimap_string_0.xml",0);
  generate_legacy_archive<
    boost::unordered_multimap<std::string,std::string>,boost::archive::xml_oarchive
  >("multimap_string_10.xml",10);
  generate_legacy_archive<
    boost::unordered_multimap<std::string,std::string>,boost::archive::xml_oarchive
  >("multimap_string_100.xml",100);

  generate_legacy_archive<
    boost::unordered_set<int>,boost::archive::xml_oarchive
  >("set_int_0.xml",0);
  generate_legacy_archive<
    boost::unordered_set<int>,boost::archive::xml_oarchive
  >("set_int_10.xml",10);
  generate_legacy_archive<
    boost::unordered_set<int>,boost::archive::xml_oarchive
  >("set_int_100.xml",100);
  generate_legacy_archive<
    boost::unordered_set<std::string>,boost::archive::xml_oarchive
  >("set_string_0.xml",0);
  generate_legacy_archive<
    boost::unordered_set<std::string>,boost::archive::xml_oarchive
  >("set_string_10.xml",10);
  generate_legacy_archive<
    boost::unordered_set<std::string>,boost::archive::xml_oarchive
  >("set_string_100.xml",100);

  generate_legacy_archive<
    boost::unordered_multiset<int>,boost::archive::xml_oarchive
  >("multiset_int_0.xml",0);
  generate_legacy_archive<
    boost::unordered_multiset<int>,boost::archive::xml_oarchive
  >("multiset_int_10.xml",10);
  generate_legacy_archive<
    boost::unordered_multiset<int>,boost::archive::xml_oarchive
  >("multiset_int_100.xml",100);
  generate_legacy_archive<
    boost::unordered_multiset<std::string>,boost::archive::xml_oarchive
  >("multiset_string_0.xml",0);
  generate_legacy_archive<
    boost::unordered_multiset<std::string>,boost::archive::xml_oarchive
  >("multiset_string_10.xml",10);
  generate_legacy_archive<
    boost::unordered_multiset<std::string>,boost::archive::xml_oarchive
  >("multiset_string_100.xml",100);

}