# Boost PolyCollection library

[![Branch](https://img.shields.io/badge/branch-master-brightgreen.svg)](https://github.com/boostorg/poly_collection/tree/master) [![CI](https://github.com/boostorg/poly_collection/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/boostorg/poly_collection/actions/workflows/ci.yml) [![Drone status](https://img.shields.io/drone/build/boostorg/poly_collection/master?server=https%3A%2F%2Fdrone.cpp.al&logo=drone&logoColor=%23CCCCCC&label=CI)](https://drone.cpp.al/boostorg/poly_collection) [![Build status](https://img.shields.io/appveyor/build/joaquintides/poly-collection/master?logo=appveyor&label=CI)](https://ci.appveyor.com/project/joaquintides/poly-collection/branch/master) [![Deps](https://img.shields.io/badge/deps-master-brightgreen.svg)](https://pdimov.github.io/boostdep-report/master/poly_collection.html)  [![Documentation](https://img.shields.io/badge/docs-master-brightgreen.svg)](https://www.boost.org/doc/libs/master/doc/html/poly_collection.html)  [![Enter the Matrix](https://img.shields.io/badge/matrix-master-brightgreen.svg)](http://www.boost.org/development/tests/master/developer/poly_collection.html)<br/>
[![Branch](https://img.shields.io/badge/branch-develop-brightgreen.svg)](https://github.com/boostorg/poly_collection/tree/develop) [![CI](https://github.com/boostorg/poly_collection/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/boostorg/poly_collection/actions/workflows/ci.yml) [![Drone status](https://img.shields.io/drone/build/boostorg/poly_collection/develop?server=https%3A%2F%2Fdrone.cpp.al&logo=drone&logoColor=%23CCCCCC&label=CI)](https://drone.cpp.al/boostorg/poly_collection) [![Build status](https://img.shields.io/appveyor/build/joaquintides/poly-collection/develop?logo=appveyor&label=CI)](https://ci.appveyor.com/project/joaquintides/poly-collection/branch/develop) [![Deps](https://img.shields.io/badge/deps-develop-brightgreen.svg)](https://pdimov.github.io/boostdep-report/develop/poly_collection.html) [![Documentation](https://img.shields.io/badge/docs-develop-brightgreen.svg)](https://www.boost.org/doc/libs/develop/doc/html/poly_collection.html) [![Enter the Matrix](https://img.shields.io/badge/matrix-develop-brightgreen.svg)](http://www.boost.org/development/tests/develop/developer/poly_collection.html)<br/>
[![BSL 1.0](https://img.shields.io/badge/license-BSL_1.0-blue.svg)](https://www.boost.org/users/license.html) <img alt="C++11 required" src="https://img.shields.io/badge/standard-C%2b%2b11-blue.svg"> <img alt="Header-only library" src="https://img.shields.io/badge/build-header--only-blue.svg">

**Boost.PolyCollection**: fast containers of polymorphic objects.

Typically, polymorphic objects cannot be stored *directly* in regular containers
and need be accessed through an indirection pointer, which introduces performance
problems related to CPU caching and branch prediction. Boost.PolyCollection
implements a
[novel data structure](http://www.boost.org/doc/html/poly_collection/an_efficient_polymorphic_data_st.html)
that is able to contiguously store polymorphic objects without such indirection,
thus providing a value-semantics user interface and better performance.
Three *polymorphic collections* are provided:

* [`boost::base_collection`](http://www.boost.org/doc/html/poly_collection/tutorial.html#poly_collection.tutorial.basics.boost_base_collection) 
* [`boost::function_collection`](http://www.boost.org/doc/html/poly_collection/tutorial.html#poly_collection.tutorial.basics.boost_function_collection)
* [`boost::any_collection`](http://www.boost.org/doc/html/poly_collection/tutorial.html#poly_collection.tutorial.basics.boost_any_collection)

dealing respectively with classic base/derived or OOP polymorphism, function wrapping
in the spirit of `std::function` and so-called
[*duck typing*](https://en.wikipedia.org/wiki/Duck_typing) as implemented by
[Boost.TypeErasure](http://www.boost.org/libs/type_erasure).

## Learn about Boost.PolyCollection

 * [Online docs](http://boost.org/libs/poly_collection)  
 * [Seminal article at bannalia.blogspot.com](http://bannalia.blogspot.com/2014/05/fast-polymorphic-collections.html)

## Install Boost.PolyCollection

* [Download Boost](https://www.boost.org/users/download/) and you're ready to go (this is a header-only library requiring no building).
* Using Conan 2: In case you don't have it yet, add an entry for Boost in your `conanfile.txt` (the example requires at least Boost 1.86):
```
[requires]
boost/[>=1.86.0]
```
<ul>If you're not using any compiled Boost library, the following will skip building altogether:</ul>

```
[options]
boost:header_only=True
```
* Using vcpkg: Execute the command
```
vcpkg install boost-poly-collection
```
* Using CMake: [Boost CMake support infrastructure](https://github.com/boostorg/cmake)
allows you to use CMake directly to download, build and consume all of Boost or
some specific libraries.

## Support

* Join the **#boost** discussion group at [cpplang.slack.com](https://cpplang.slack.com/)
([ask for an invite](https://cppalliance.org/slack/) if you’re not a member of this workspace yet)
* Ask in the [Boost Users mailing list](https://lists.boost.org/mailman/listinfo.cgi/boost-users)
(add the `[poly_collection]` tag at the beginning of the subject line)
* [File an issue](https://github.com/boostorg/poly_collection/issues)

## Contribute

* [Pull requests](https://github.com/boostorg/poly_collection/pulls) against **develop** branch are most welcome.
Note that by submitting patches you agree to license your modifications under the [Boost Software License, Version 1.0](http://www.boost.org/LICENSE_1_0.txt).
