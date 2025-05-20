+++
title = "Using a Result"
description = "Using a C Result"
weight = 30
+++

This models [the earlier C++ example of use]({{% relref "/experimental/worked-example/implicit-construction" %}}),
and its C equivalent isn't much more verbose thanks to our helper typedefs and macros:

{{% snippet "c_api2.cpp" "using" %}}

For this to link, the `BOOST_OUTCOME_C_DECLARE_RESULT_SYSTEM_FROM_ENUM` macro needs to be
compiled at least once within C++ within the final binary to emit the extern
functions needed by C.
