<a href="http://mc-stan.org">
<img src="https://github.com/stan-dev/stan/blob/master/logos/stanlogo-main.png?raw=true" alt="Stan Logo"/>
</a>


The <b>Stan Math Library</b> is a C++, reverse-mode automatic differentiation library designed to be usable, extensive and extensible, efficient, scalable, stable, portable, and redistributable in order to facilitate the construction and utilization of algorithms that utilize derivatives.

Home Page
---------
Stan's home page, with links to everything you'll need to use the Stan Math Library is:

[http://mc-stan.org/](http://mc-stan.org/)

Licensing
---------
The Stan Math Library is licensed under new BSD.


Installation
------------
The Stan Math Library is a header-only C++ library. To use the Stan Math Library, add the current directory to the source path (in `g++` and `clang++` it is the `-I <math directory>` option) and add this to your source file: `#include <stan/math.hpp>`.

Running Tests
-------------
Running tests in the Stan Math Library require these libraries:

- Boost
- Eigen
- Google Test
- CppLint (optional)

If the Stan Math Library has been cloned as part of the Stan library, no configuration is necessary.

If the Stan Math Library is outside of the Stan library, the following variables have to be provided to make in either `make/local` or `~/.config/stan/make.local`:

- `STANAPI_HOME` (the home directory of the Stan library)

or

- `EIGEN`
- `BOOST`
- `GTEST`
- `CPPLINT` (optional)

If you see this message:

```
------------------------------------------------------------
make generate-tests -s
test/prob/generate_tests.cpp:9:10: fatal error: 'boost/algorithm/string.hpp' file not found
#include <boost/algorithm/string.hpp>
         ^
1 error generated.
make: *** [test/prob/generate_tests] Error 1
make generate-tests -s failed
```

the library paths have not been configured correctly.
