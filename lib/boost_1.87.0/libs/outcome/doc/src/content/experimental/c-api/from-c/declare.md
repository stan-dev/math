+++
title = "Declare a Result"
description = "Declaring a C Result"
weight = 20
+++

{{% snippet "c_api2.cpp" "preamble" %}}

The key to making C programming easy is to alias the long complex things
into short easy thing. Obviously `SUCCESS(expr)` and `FAILURE(expr)` is too
generic, but for the purposes of this documentation it makes thing easier.
