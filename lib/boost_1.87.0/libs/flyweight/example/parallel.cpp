/* Boost.Flyweight example of parallel tokenization.
 *
 * Copyright 2024 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See http://www.boost.org/libs/flyweight for library home page.
 */

#include <boost/flyweight.hpp>
#include <boost/flyweight/concurrent_factory.hpp>
#include <boost/flyweight/no_locking.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

/* Handcrafted tokenizer for sequences of alphabetic characters */

inline bool match(char ch)
{
 return (ch>='a' && ch<='z') || (ch>='A' && ch<='Z');
}

template<typename ForwardIterator,typename F>
void tokenize(ForwardIterator first,ForwardIterator last,F f)
{
  goto start;
  for(;;)
  {
    for(;;){
      ++first;

    start:
      if(first==last)return;
      if(match(*first))break;
    }

    auto begin_word=first;
    for(;;){
      if(++first==last||!match(*first)){
        f(begin_word,first);
        if(first==last)return;
        else           break;
      }
    }
  }
}

/* Tokenize a string into words in parallel and store the results into a
 * std::vector<String>, String being std::string or a flyweight type.
 */

template<typename String>
void parse(const std::string& in,const char* type_name,std::size_t num_threads)
{
  using namespace std::chrono;
  using string_iterator=std::string::const_iterator;

  auto t1=steady_clock::now();

  /* Divide input in num_threads chunks, taking care that boundaries are not
   * placed in the middle of a token.
   */

  std::vector<string_iterator> boundaries(num_threads+1);

  boundaries[0]=in.begin();
  for(std::size_t i=0;i<num_threads;++i){
    auto& it=boundaries[i+1];
    it=boundaries[i]+(in.end()-boundaries[i])/(num_threads-i);
    while(it!=in.end()&&match(*it))++it;
  }

  /* do a first pass to precalculate # of words produced by each thread */

  std::vector<std::thread> threads(num_threads);
  std::vector<std::size_t> partial_num_words(num_threads);

  for(std::size_t i=0;i<num_threads;++i){
    threads[i]=std::thread([&,i]{
      std::size_t s=0;
      tokenize(
        boundaries[i],boundaries[i+1],
        [&](string_iterator,string_iterator){++s;});
      partial_num_words[i]=s;
    });
  }

  std::size_t              num_words=0;
  std::vector<std::size_t> thread_output_starts(num_threads);

  for(std::size_t i=0;i<num_threads;++i){
    threads[i].join();
    thread_output_starts[i]=num_words;
    num_words+=partial_num_words[i];
  }

  /* do a second pass, this time populating the result vector */

  std::vector<String> words(num_words,String());

  for(std::size_t i=0;i<num_threads;++i){
    threads[i]=std::thread([&,i]{
      auto out=words.begin()+thread_output_starts[i];
      tokenize(
        boundaries[i],boundaries[i+1],
        [&](string_iterator first,string_iterator last){
          *out++=String(first,last);
        });
    });
  }
  for(std::size_t i=0;i<num_threads;++i){threads[i].join();}

  auto t2=steady_clock::now();

  std::cout
    <<std::setw(20)<<type_name<<", "<<num_threads<<" thread(s): "
    <<num_words<<" words, "
    <<std::setw(9)<<duration_cast<duration<double>>(t2-t1).count()<< " s\n";
}

/* accept a file and parse it with std::string and various flyweight types */

int main(int argc, char** argv)
{
  using namespace boost::flyweights;
  using regular_flyweight=flyweight<std::string>;
  using concurrent_flyweight=flyweight<
    std::string,
    concurrent_factory<>,
    no_locking,
    no_tracking
  >;

  if(argc<2){
    std::cout<<"specify a file\n";
    std::exit(EXIT_FAILURE);
  }

  std::ifstream is(argv[1]);
  if(!is)
  {
    std::cout<<"can't open "<<argv[1]<<"\n";
    std::exit(EXIT_FAILURE);
  }

  std::string in(
    std::istreambuf_iterator<char>(is),std::istreambuf_iterator<char>{});

  parse<std::string>(in,"std::string",1);
  parse<std::string>(in,"std::string",8);
  parse<regular_flyweight>(in,"regular flyweight",1);
  parse<regular_flyweight>(in,"regular flyweight",8);
  parse<concurrent_flyweight>(in,"concurrent flyweight",1);
  parse<concurrent_flyweight>(in,"concurrent flyweight",8);
}
