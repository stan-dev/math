/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------*/

#include <iostream>
#include <string>

#define SUNDIALS_HOST_DEVICE
#define SUNDIALS_DEVICE_INLINE
#include "sundials/sundials_reductions.hpp"

using namespace sundials::reductions;
using namespace sundials::reductions::impl;

int testPlusWithInts()
{
  const std::string testStr = "Running testPlusWithInts";

  std::cout << testStr;

  int a{5};
  int b{1};
  int c{0};
  c = plus<int>{}(a,b);

  bool pass = c == a + b;

  if (pass)
    std::cout << " -- passed\n";
  else
    std::cout << " -- FAILED\n";

  return !pass;
}

int testPlusWithDoubles()
{
  const std::string testStr = "Running testPlusWithDoubles";

  std::cout << testStr;

  double a{5.0};
  double b{1.0};
  double c{0.0};
  c = plus<double>{}(a,b);

  bool pass = c == a + b;

  if (pass)
    std::cout << " -- passed\n";
  else
    std::cout << " -- FAILED\n";

  return !pass;
}

int testMaximumWithInts()
{
  const std::string testStr = "Running testMaximumWithInts";

  std::cout << testStr;

  int a{5};
  int b{1};
  int c{0};
  c = maximum<int>{}(a,b);
  bool pass1 = c == a;

  a = 0;
  c = maximum<int>{}(a,b);
  bool pass2 = c == b;

  bool pass = pass1 && pass2;
  if (pass)
    std::cout << " -- passed\n";
  else
    std::cout << " -- FAILED\n";

  return !pass;
}

int testMaximumWithDoubles()
{
  const std::string testStr = "Running testMaximumWithDoubles";

  std::cout << testStr;

  double a{5.0};
  double b{1.0};
  double c{0.0};
  c = maximum<double>{}(a,b);
  bool pass1 = c == a;

  a = 0.0;
  c = maximum<double>{}(a,b);
  bool pass2 = c == b;

  bool pass = pass1 && pass2;
  if (pass)
    std::cout << " -- passed\n";
  else
    std::cout << " -- FAILED\n";

  return !pass;
}

int testMinimumWithInts()
{
  const std::string testStr = "Running testMinimumWithInts";

  std::cout << testStr;

  int a{5};
  int b{1};
  int c{0};
  c = minimum<int>{}(a,b);
  bool pass1 = c == b;

  a = 0;
  c = minimum<int>{}(a,b);
  bool pass2 = c == a;

  bool pass = pass1 && pass2;
  if (pass)
    std::cout << " -- passed\n";
  else
    std::cout << " -- FAILED\n";

  return !pass;
}

int testMinimumWithDoubles()
{
  const std::string testStr = "Running testMinimumWithDoubles";

  std::cout << testStr;

  double a{5.0};
  double b{1.0};
  double c{0.0};
  c = minimum<double>{}(a,b);
  bool pass1 = c == b;

  a = 0.0;
  c = minimum<double>{}(a,b);
  bool pass2 = c == a;

  bool pass = pass1 && pass2;
  if (pass)
    std::cout << " -- passed\n";
  else
    std::cout << " -- FAILED\n";

  return !pass;
}

int main()
{
  int fails{0};

  std::cout << "Testing sundials::reductions operations\n";

  fails += testPlusWithInts();
  fails += testPlusWithDoubles();
  fails += testMaximumWithInts();
  fails += testMaximumWithDoubles();
  fails += testMinimumWithInts();
  fails += testMinimumWithDoubles();

  if (fails)
    std::cout << "FAIL: " << fails << " tests failed\n";
  else
    std::cout << "SUCCESS: all tests passed\n";

  return fails;
}

