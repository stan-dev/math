//  Copyright (c) 2003 Vladimir Prus
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)
//
//  http://www.boost.org
//

#include <a.h>

int main() {
#ifdef RELEASE
  release();
#else
  debug();
#endif
  return 0;
}
