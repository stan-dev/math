+++
title = "TRY a C Result"
description = "Operation TRY on a C Result"
weight = 40
+++

Thanks to much of the magic of {{< api "BOOST_OUTCOME_TRY(var, expr)" >}} being implemented
using C preprocessor metaprogramming, we can offer a very similar experience for the
C try operation and without needing anything compiled in C++ as support functions:

{{% snippet "c_api2.cpp" "try" %}}

The principle difference is that you can specify a cleanup routine to perform if
failure is encountered. This is especially useful in C, which has no stack unwinding.

Also due to lack of type sugaring and user defined implicit conversions, if your
callers result type isn't your callee's, you may need to specify what your caller's
result type is so the error state can be correctly propagated.
