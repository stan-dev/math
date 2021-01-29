# ![Boost.Geometry](doc/other/logo/logo_bkg.png)

Boost.Geometry, part of collection of the [Boost C++ Libraries](http://github.com/boostorg), defines concepts, primitives and algorithms for solving geometry problems.

[![Licence](https://img.shields.io/badge/license-boost-4480cc.png)](http://www.boost.org/LICENSE_1_0.txt)
[![Documentation](https://img.shields.io/badge/-documentation-4480cc.png)](http://boost.org/libs/geometry)
[![Wiki](https://img.shields.io/badge/-wiki-4480cc.png)](https://github.com/boostorg/geometry/wiki)
[![Mailing List](https://img.shields.io/badge/-mailing%20list-4eb899.png)](http://lists.boost.org/geometry/)
[![Chat](https://badges.gitter.im/boostorg/geometry.png)](https://gitter.im/boostorg/geometry?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

| :warning: **CAUTION**: Boost.Geometry in Boost 1.73 deprecates support for the C++03 and will require C++14 from Boost 1.75 onwards (see issue [#590](https://github.com/boostorg/geometry/issues/590)). |
| --- |

### Test results

 Branch     | Build         | Coverage       | Regression | Documentation
------------|---------------|----------------|------------|--------------
**develop** | [![status](https://circleci.com/gh/boostorg/geometry/tree/develop.svg?style=shield)](https://circleci.com/gh/boostorg/geometry/tree/develop) <br> [![status](https://github.com/boostorg/geometry/workflows/clang-test-minimal/badge.svg?branch=develop)](https://github.com/boostorg/geometry/actions?query=branch:develop+workflow:clang-test-minimal) <br> [![status](https://github.com/boostorg/geometry/workflows/gcc-test-minimal/badge.svg?branch=develop)](https://github.com/boostorg/geometry/actions?query=branch:develop+workflow:gcc-test-minimal) <br> [![status](https://github.com/boostorg/geometry/workflows/msvc-test-minimal/badge.svg?branch=develop)](https://github.com/boostorg/geometry/actions?query=branch:develop+workflow:msvc-test-minimal) | [![status](https://coveralls.io/repos/github/boostorg/geometry/badge.svg?branch=develop)](https://coveralls.io/github/boostorg/geometry?branch=develop) | [![geometry](https://img.shields.io/badge/-geometry-4480cc.png)](http://www.boost.org/development/tests/develop/developer/geometry.html) [![index](https://img.shields.io/badge/-index-4480cc.png)](http://www.boost.org/development/tests/develop/developer/geometry-index.html) [![extensions](https://img.shields.io/badge/-extensions-4480cc.png)](http://www.boost.org/development/tests/develop/developer/geometry-extensions.html) | [![](https://github.com/boostorg/geometry/workflows/documentation/badge.svg?branch=develop)](https://github.com/boostorg/geometry/actions?query=branch:develop+workflow:documentation)
**master**  | [![status](https://circleci.com/gh/boostorg/geometry/tree/master.svg?style=shield)](https://circleci.com/gh/boostorg/geometry/tree/master)  <br> [![status](https://github.com/boostorg/geometry/workflows/clang-test-minimal/badge.svg?branch=master)](https://github.com/boostorg/geometry/actions?query=branch:master+workflow:clang-test-minimal) <br> [![status](https://github.com/boostorg/geometry/workflows/gcc-test-minimal/badge.svg?branch=master)](https://github.com/boostorg/geometry/actions?query=branch:master+workflow:gcc-test-minimal) <br> [![status](https://github.com/boostorg/geometry/workflows/msvc-test-minimal/badge.svg?branch=master)](https://github.com/boostorg/geometry/actions?query=branch:master+workflow:msvc-test-minimal) | [![status](https://coveralls.io/repos/github/boostorg/geometry/badge.svg?branch=master)](https://coveralls.io/github/boostorg/geometry?branch=master)   | [![geometry](https://img.shields.io/badge/-geometry-4480cc.png)](http://www.boost.org/development/tests/master/developer/geometry.html) [![index](https://img.shields.io/badge/-index-4480cc.png)](http://www.boost.org/development/tests/master/developer/geometry-index.html) | [![](https://github.com/boostorg/geometry/workflows/documentation/badge.svg?branch=master)](https://github.com/boostorg/geometry/actions?query=branch:master+workflow:documentation)

### Directories

* **doc** - QuickBook documentation sources
* **example** - Boost.Geometry examples
* **_extensions_** - examples and tests for the extensions - _develop branch_
* **include** - the sourcecode of Boost.Geometry
* **index** - examples and tests for the Spatial Index
* **meta** - library metadata
* **test** - Boost.Geometry unit tests
