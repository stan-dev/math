+++
title = "`static void on_result_move_construction(T *, U &&) noexcept`"
description = "Hook invoked by the converting move constructors of `basic_result`."
categories = ["observer-policies"]
weight = 450
+++

One of the constructor hooks for {{% api "basic_result<T, E, NoValuePolicy>" %}}, generally invoked by the converting move constructors of `basic_result` (NOT the standard move constructor). See each constructor's documentation to see which specific hook it invokes.

*Requires*: Always available.

*Guarantees*: Never throws an exception.
