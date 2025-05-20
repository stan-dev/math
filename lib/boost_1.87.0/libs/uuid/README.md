# Boost.Uuid

Boost.Uuid, part of [Boost C++ Libraries](http://boost.org),
provides a C++ implementation of Universally Unique Identifiers (UUID) as described in [RFC 4122](https://datatracker.ietf.org/doc/rfc4122/) and [RFC 9562](https://datatracker.ietf.org/doc/rfc9562/).

See [the documentation](http://boost.org/libs/uuid) for more information.

## License

Distributed under the [Boost Software License, Version 1.0](https://www.boost.org/LICENSE_1_0.txt).

## Properties

* C++11 (since Boost 1.86.0)
* Header-only

## Current Status

Branch          | Github Actions | Appveyor | Dependencies | Documentation | Test Matrix |
:-------------: | -------------- | -------- | ------------ | ------------- | ----------- |
[`master`](https://github.com/boostorg/uuid/tree/master) | [![Build Status](https://github.com/boostorg/uuid/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/boostorg/uuid/actions?query=branch:master) | [![Build status](https://ci.appveyor.com/api/projects/status/rmp9xmse2b6urkjv/branch/master?svg=true)](https://ci.appveyor.com/project/cppalliance/uuid/branch/master) | [![Deps](https://img.shields.io/badge/deps-master-brightgreen.svg)](https://pdimov.github.io/boostdep-report/master/uuid.html) | [![Documentation](https://img.shields.io/badge/docs-master-brightgreen.svg)](http://www.boost.org/doc/libs/master/doc/html/uuid.html) | [![Enter the Matrix](https://img.shields.io/badge/matrix-master-brightgreen.svg)](http://www.boost.org/development/tests/master/developer/uuid.html)
[`develop`](https://github.com/boostorg/uuid/tree/develop) | [![Build Status](https://github.com/boostorg/uuid/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/boostorg/uuid/actions?query=branch:develop) | [![Build status](https://ci.appveyor.com/api/projects/status/rmp9xmse2b6urkjv/branch/develop?svg=true)](https://ci.appveyor.com/project/cppalliance/uuid/branch/develop) | [![Deps](https://img.shields.io/badge/deps-develop-brightgreen.svg)](https://pdimov.github.io/boostdep-report/develop/uuid.html) | [![Documentation](https://img.shields.io/badge/docs-develop-brightgreen.svg)](http://www.boost.org/doc/libs/develop/doc/html/uuid.html) | [![Enter the Matrix](https://img.shields.io/badge/matrix-develop-brightgreen.svg)](http://www.boost.org/development/tests/develop/developer/uuid.html)

## More Information

* [Ask questions](http://stackoverflow.com/questions/ask?tags=c%2B%2B,boost,boost-uuid)
* [Report bugs](https://github.com/boostorg/uuid/issues): Be sure to mention Boost version, platform and compiler you're using. A small compilable code sample to reproduce the problem is always good as well.
* Submit your patches as pull requests against the **develop** branch. Note that by submitting patches you agree to license your modifications under the [Boost Software License, Version 1.0](https://www.boost.org/LICENSE_1_0.txt).
* Discussions about the library are held on the [Boost developers mailing list](http://www.boost.org/community/groups.html#main). Be sure to read the [discussion policy](http://www.boost.org/community/policy.html) before posting and add the `[uuid]` tag at the beginning of the subject line.

## Code Example - UUID Generation

```cpp
//  mkuuid.cpp example

#include <boost/uuid.hpp>
#include <iostream>

int main()
{
    boost::uuids::random_generator gen;
    std::cout << gen() << std::endl;
}
```

```shell
$ clang++ -Wall -Wextra -std=c++11 -O2 mkuuid.cpp -o mkuuid
$ ./mkuuid
2c186eb0-89cf-4a3c-9b97-86db1670d5f4
$ ./mkuuid
a9d3fbb9-0383-4389-a8a8-61f6629f90b6
```
