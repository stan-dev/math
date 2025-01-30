# Boost.Flyweight

[![Branch](https://img.shields.io/badge/branch-master-brightgreen.svg)](https://github.com/boostorg/flyweight/tree/master) [![CI](https://github.com/boostorg/flyweight/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/boostorg/flyweight/actions/workflows/ci.yml) [![Drone status](https://img.shields.io/drone/build/boostorg/flyweight/master?server=https%3A%2F%2Fdrone.cpp.al&logo=drone&logoColor=%23CCCCCC&label=CI)](https://drone.cpp.al/boostorg/flyweight)  [![Deps](https://img.shields.io/badge/deps-master-brightgreen.svg)](https://pdimov.github.io/boostdep-report/master/flyweight.html)  [![Documentation](https://img.shields.io/badge/docs-master-brightgreen.svg)](https://www.boost.org/doc/libs/master/libs/flyweight)  [![Enter the Matrix](https://img.shields.io/badge/matrix-master-brightgreen.svg)](http://www.boost.org/development/tests/master/developer/flyweight.html)<br/>
[![Branch](https://img.shields.io/badge/branch-develop-brightgreen.svg)](https://github.com/boostorg/flyweight/tree/develop) [![CI](https://github.com/boostorg/flyweight/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/boostorg/flyweight/actions/workflows/ci.yml) [![Drone status](https://img.shields.io/drone/build/boostorg/flyweight/develop?server=https%3A%2F%2Fdrone.cpp.al&logo=drone&logoColor=%23CCCCCC&label=CI)](https://drone.cpp.al/boostorg/flyweight)  [![Deps](https://img.shields.io/badge/deps-develop-brightgreen.svg)](https://pdimov.github.io/boostdep-report/develop/flyweight.html) [![Documentation](https://img.shields.io/badge/docs-develop-brightgreen.svg)](https://www.boost.org/doc/libs/develop/libs/flyweight) [![Enter the Matrix](https://img.shields.io/badge/matrix-develop-brightgreen.svg)](http://www.boost.org/development/tests/develop/developer/flyweight.html)<br/>
[![BSL 1.0](https://img.shields.io/badge/license-BSL_1.0-blue.svg)](https://www.boost.org/users/license.html) <img alt="Header-only library" src="https://img.shields.io/badge/build-header--only-blue.svg">

Flyweights are small-sized handle classes granting constant access to shared common data, thus allowing for the management
of large amounts of entities within reasonable memory limits. Boost.Flyweight makes it easy to use this common programming
idiom by providing the class template `flyweight<T>`, which acts as a drop-in replacement for `const T`.

## Learn about Boost.Flyweight

* [Online documentation](https://boost.org/libs/flyweight)

## Install Boost.Flyweight

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
vcpkg install boost-flyweight
```
* Using CMake: [Boost CMake support infrastructure](https://github.com/boostorg/cmake)
allows you to use CMake directly to download, build and consume all of Boost or
some specific libraries.

## Support

* Join the **#boost** discussion group at [cpplang.slack.com](https://cpplang.slack.com/)
([ask for an invite](https://cppalliance.org/slack/) if you’re not a member of this workspace yet)
* Ask in the [Boost Users mailing list](https://lists.boost.org/mailman/listinfo.cgi/boost-users)
(add the `[flyweight]` tag at the beginning of the subject line)
* [File an issue](https://github.com/boostorg/flyweight/issues)

## Contribute

* [Pull requests](https://github.com/boostorg/flyweight/pulls) against **develop** branch are most welcome.
Note that by submitting patches you agree to license your modifications under the [Boost Software License, Version 1.0](http://www.boost.org/LICENSE_1_0.txt).
