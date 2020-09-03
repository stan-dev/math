## Windows Development Tips {#windows_tips}

It seems to be challenging these days to get a working Unix compilation toolchain on Windows, moreso than it was in the past. Here are some notes on what worked for one developer.

The complilation shell environment is [MSYS2 MinGW 32-bit](http://www.msys2.org/). However, any git commands are executed separately in [git for Windows](https://git-scm.com/download/win)' git bash. All commands herein have been run from the top-most `math` folder of this repo.

```
$ which g++
/c/Users/<user name>/AppData/Local/mingw32/bin/g++

$ g++ --version
g++.exe (i686-posix-dwarf-rev0, Built by MinGW-W64 project) 7.2.0
Copyright (C) 2017 Free Software Foundation, Inc.
```

This compiler came from the [MinGW-W64 builds](https://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/) project. To get the compiler to invoke with just `g++`, one may have to [edit the environment variables for the local account](https://superuser.com/a/989665/561353) to append `%USERPROFILE%\AppData\Local\mingw32\bin` to the `PATH`, set `MSYS2_PATH_TYPE=inherit`, and follow with a reboot.

The make and python based build structure referenced in the [stan developer overview](https://github.com/stan-dev/stan/wiki/Developer-process-overview) is replicated in the stan-math repo, so one can use the instructions there to get an idea of what commands the build system supports, like `make test-headers` or `runTests.py`. Be aware that these commands can take a long time, and that stan-math's automated build system, Travis, is running older compilers like GCC 4.9. Therefore, new compiler bugs may be encountered deep in the middle of the tests.

```
$ which python
/mingw32/bin/python

$ pacman -Ss mingw-w64-i686-python2 | grep installed
mingw32/mingw-w64-i686-python2 2.7.12-1 [installed]

$ python --version
Python 2.7.12
```
The shebang at the beginning of `runTests.py` may need to be shortened to just say `python`, rather than a full path.
```
#!/usr/bin/python
... rest of listing omitted ...
```
The build system supports making and running individual tests.
```
$ make test/unit/math/rev/core/var_test.exe
make: 'test/unit/math/rev/core/var_test.exe' is up to date.

$ ./test/unit/math/rev/core/var_test.exe
Running main() from gtest_main.cc
... lines omitted ...
[==========] 16 tests from 1 test case ran. (0 ms total)
[  PASSED  ] 16 tests.
```
Autodetection in the build system may fail, so one may need to create a `local` file under `make` to help. The wiki authors do recognize that technically `CC` should be `gcc`.
```
CC=g++
CXX=g++
BIT=32
```
If the test suite isn't built first, client code using numerical integration routines such as `integrate_ode_bdf` may fail to link because the libraries haven't been built yet. Obviously these files and directories would need to be added to the downstream project makefile's LDLIBS variable, or equivalent.
```
$ make lib/sundials_3.3.0/lib/libsundials_cvodes.a lib/sundials_3.3.0/lib/libsundials_nvecserial.a
make: 'lib/sundials_3.3.0/lib/libsundials_cvodes.a' is up to date.
make: 'lib/sundials_3.3.0/lib/libsundials_nvecserial.a' is up to date.
```
