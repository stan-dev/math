//          Copyright Alain Miniussi 2023
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

// Authors: Alain Miniussi

#ifndef BOOST_MPI_MPI_TEST_UTILS_HPP
#define BOOST_MPI_MPI_TEST_UTILS_HPP

#include <iostream>

inline
void
check_failed(bool cond, std::string msg, int& failed) {
  if (!cond) {
    std::cerr << "FAILED: " << msg << '\n';
    ++failed;
  } else {
    std::cerr << "PASSED: " << msg << '\n';
  }
}

template<typename T>
void
check_failed(T cond, std::string msg, int& failed) {
  check_failed(bool(cond), msg, failed);
}

inline
void
count_failed(int nfailed, std::string msg, int& failed) {
  if (nfailed > 0) {
    std::cerr << "FAILED: " << nfailed << " in " << msg << '\n';
    failed += nfailed;
  } else {
    std::cerr << "PASSED: " << msg << '\n';
  }
}

template<class T>
void
count_failed(T nfailed, std::string msg, int& failed) {
  count_failed(int(nfailed), msg, failed);
}

#define BOOST_MPI_CHECK(cond, failed)        check_failed(cond, #cond, failed);
#define BOOST_MPI_COUNT_FAILED(fct,  failed) count_failed(fct, #fct, failed);

#endif
