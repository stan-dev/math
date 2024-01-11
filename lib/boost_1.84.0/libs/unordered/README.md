# Boost.Unordered

[![Branch](https://img.shields.io/badge/branch-master-brightgreen.svg)](https://github.com/boostorg/unordered/tree/master) [![CI](https://github.com/boostorg/unordered/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/boostorg/unordered/actions/workflows/ci.yml) [![Drone status](https://img.shields.io/drone/build/boostorg/unordered/master?server=https%3A%2F%2Fdrone.cpp.al&logo=drone&logoColor=%23CCCCCC&label=CI)](https://drone.cpp.al/boostorg/unordered) [![Build status](https://img.shields.io/appveyor/build/cppalliance/unordered/master?logo=appveyor&label=CI)](https://ci.appveyor.com/project/cppalliance/unordered/branch/master)  [![codecov](https://codecov.io/gh/boostorg/unordered/branch/master/graph/badge.svg)](https://codecov.io/gh/boostorg/unordered/branch/master)  [![Deps](https://img.shields.io/badge/deps-master-brightgreen.svg)](https://pdimov.github.io/boostdep-report/master/unordered.html)  [![Documentation](https://img.shields.io/badge/docs-master-brightgreen.svg)](https://www.boost.org/doc/libs/master/libs/unordered/doc/html/unordered.html)  [![Enter the Matrix](https://img.shields.io/badge/matrix-master-brightgreen.svg)](http://www.boost.org/development/tests/master/developer/unordered.html)<br/>
[![Branch](https://img.shields.io/badge/branch-develop-brightgreen.svg)](https://github.com/boostorg/unordered/tree/develop) [![CI](https://github.com/boostorg/unordered/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/boostorg/unordered/actions/workflows/ci.yml) [![Drone status](https://img.shields.io/drone/build/boostorg/unordered/develop?server=https%3A%2F%2Fdrone.cpp.al&logo=drone&logoColor=%23CCCCCC&label=CI)](https://drone.cpp.al/boostorg/unordered) [![Build status](https://img.shields.io/appveyor/build/cppalliance/unordered/master?logo=appveyor&label=CI)](https://ci.appveyor.com/project/cppalliance/unordered/branch/develop) [![codecov](https://codecov.io/gh/boostorg/unordered/branch/develop/graph/badge.svg)](https://codecov.io/gh/boostorg/unordered/branch/develop) [![Deps](https://img.shields.io/badge/deps-develop-brightgreen.svg)](https://pdimov.github.io/boostdep-report/develop/unordered.html) [![Documentation](https://img.shields.io/badge/docs-develop-brightgreen.svg)](https://www.boost.org/doc/libs/develop/libs/unordered/doc/html/unordered.html) [![Enter the Matrix](https://img.shields.io/badge/matrix-develop-brightgreen.svg)](http://www.boost.org/development/tests/develop/developer/unordered.html)<br/>
[![BSL 1.0](https://img.shields.io/badge/license-BSL_1.0-blue.svg)](https://www.boost.org/users/license.html) <img alt="C++11 required" src="https://img.shields.io/badge/standard-C%2b%2b11-blue.svg"> <img alt="Header-only library" src="https://img.shields.io/badge/build-header--only-blue.svg">

Boost.Unordered offers a catalog of hash containers with different standards compliance levels, performances and intented usage scenarios:

**`boost::unordered_set` `boost::unordered_map` `boost::unordered_multiset` `boost::unordered_multimap`**

<ul>Fully conformant implementations of <code>std::unordered_[multi](set|map)</code>,
but faster and up to the latest revisions of the standard even if you're working in an older version of C++ (heterogeneous lookup,
<code>try_emplace</code>, <code>contains</code>, etc.)</ul>

**`boost::unordered_flat_set` `boost::unordered_flat_map`**

<ul>The fastest of the lot. Based on open addressing, these containers slightly
deviate from the standard in exchange for top performance.</ul>

**`boost::unordered_node_set` `boost::unordered_node_map`**

<ul>Variations of <code>boost::unordered_flat_(set|map)</code> providing pointer stability.</ul>

**`boost::concurrent_flat_set` `boost::concurrent_flat_map`**

<ul>High performance for multithreaded scenarios. Introducing a new non-standard, iterator-free API.</ul>

## Learn about Boost.Unordered

* [Online documentation](https://boost.org/libs/unordered)
* [Some benchmarks](https://github.com/boostorg/boost_unordered_benchmarks)
* Technical articles on Boost.Unordered internal design:
  * [Advancing the state of the art for `std::unordered_map` implementations](https://bannalia.blogspot.com/2022/06/advancing-state-of-art-for.html)
  * [Inside `boost::unordered_flat_map`](https://bannalia.blogspot.com/2022/11/inside-boostunorderedflatmap.html)
  * [Inside `boost::concurrent_flat_map`](https://bannalia.blogspot.com/2023/07/inside-boostconcurrentflatmap.html)
  * [Bulk visitation in `boost::concurrent_flat_map`](https://bannalia.blogspot.com/2023/10/bulk-visitation-in-boostconcurrentflatm.html)

## Get the library

Boost.Unordered can be installed in a number of ways:

* [Download Boost](https://www.boost.org/users/download/) and you're ready to go (this is a header-only library requiring no building).
* Using Conan 2: In case you don't have it yet, add an entry for Boost in your `conanfile.txt` (the example requires at least Boost 1.83):
```
[requires]
boost/[>=1.83.0]
```
<ul>If you're not using any compiled Boost library, the following will skip building altogether:</ul>

```
[options]
boost:header_only=True
```
* Using vcpkg: Execute the command
```
vcpkg install boost-unordered
```
* Using CMake: [Boost CMake support infrastructure](https://github.com/boostorg/cmake)
allows you to use CMake directly to download, build and consume all of Boost or
some specific libraries.

## Support

* Join the **#boost-unordered** discussion group at [cpplang.slack.com](https://cpplang.slack.com/)
([ask for an invite](https://cppalliance.org/slack/) if you’re not a member of this workspace yet)
* Ask in the [Boost Users mailing list](https://lists.boost.org/mailman/listinfo.cgi/boost-users)
(add the `[unordered]` tag at the beginning of the subject line)
* [File an issue](https://github.com/boostorg/unordered/issues)

## Contribute

* [Pull requests](https://github.com/boostorg/unordered/pulls) against **develop** branch are most welcome.
Note that by submitting patches you agree to license your modifications under the [Boost Software License, Version 1.0](http://www.boost.org/LICENSE_1_0.txt).
