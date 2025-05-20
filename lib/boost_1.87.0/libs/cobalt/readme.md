# boost.cobalt

This library provides a set of easy to use coroutine primitives & utilities running on top of boost.asio.
These will be of interest for applications that perform a lot of IO that want to not block unnecessarily,
yet still want to have linear & readable code (i..e. avoid callbacks).

A minimum of Boost 1.82 is necessary as the ASIO in that version has needed support. C++ 20 is needed for C++ coroutines.

Below is a showcase of features, if you're new to coroutines or asynchronous programming, please see the [primer](https://www.boost.org/doc/libs/master/libs/cobalt/doc/html/index.html#coroutine_primer).

The assumptions are:

 - `io_context` is the execution_context of choice.
 - If `asio::io_context` is the executor, no more than one kernel thread executes within it at a time.
 - Eager execution is the way to go.
 - A thread created with promise is only using promise stuff.

## Entry points

```cpp
// a single threaded main running on an io_context
cobalt::main co_main(int argc, char ** argv)
{
    // wrapper around asio::steady_timer
    asio::steady_timer tim{co_await cobalt::this_coro::executor};
    tim.expires_after(std::chrono::milliseconds(100));

    co_await tim.async_wait(cobalt::use_op);
    co_return 0;
}
```

That is, [`main`](doc/reference/main.adoc) runs on a single threaded `io_context`.

It also hooks up signals, so that things like `Ctrl+C` get forwarded as cancellations automatically

Alternatively, [`run`](doc/reference/run.adoc) can be used manually.

```cpp
cobalt::task<int> main_func()
{
    asio::steady_timer tim{co_await cobalt::this_coro::executor};
    tim.expires_after(std::chrono::milliseconds(100));

    co_await tim.async_wait(cobalt::use_op);
    co_return 0;
}


int main(int argc, char ** argv)
{
    return run(main_func());
}
```

## Promises

The core primitive for creating your own functions is [`cobalt::promise<T>`](doc/reference/promise.adoc).
It is eager, i.e. it starts execution immediately, before you `co_await`.

```cpp
cobalt::promise<void> test()
{
    printf("test-1\n");
    asio::steady_timer tim{co_await cobalt::this_coro::executor};
    tim.expires_after(std::chrono::milliseconds(100));
    co_await tim.async_wait(cobalt::use_op);
    printf("test-2\n");
}

cobalt::main co_main(int argc, char ** argv)
{
    printf("main-1\n");
    auto tt = test();
    printf("main-2\n");
    co_await tt;
    printf("main-3\n");
    return 0;
}
```

The output of the above will be:

```cpp
main-1
test-1
main-2
test-2
main-3
```

Unlike ops, returned by .wait, the promise can be disregarded; disregarding the promise does not cancel it, but rather detaches is. This makes it easy to 
spin up multiple tasks to run in parallel. In order to avoid accidental detaching the promise type uses `nodiscard` unless one uses `+` to detach it:

```cpp
cobalt::promise<void> my_task();

cobalt::main co_main()
{
    // warns & cancels the task
    my_task();
    // ok
    +my_task();
    co_return 0;
}
```

## Task

A [`task`](doc/reference/task.adoc) is a lazy alternative to a promise, that can be spawned onto or `co_await`ed on another executor.

An `cobalt::task` can also be used with `spawn` to turn it into an asio operation.

## Generator

A [`generator`](doc/reference/generator.adoc) is a coroutine that produces a series of values instead of one, but otherwise similar to `promise`.

```cpp
cobalt::generator<int> test()
{
  printf("test-1\n");
  co_yield 1;
  printf("test-2\n");
  co_yield 2;
  printf("test-3\n");
  co_return 3;
}

cobalt::main co_main(int argc, char ** argv)
{
    printf("main-1\n");
    auto tt = test();
    printf("main-2\n");
    i = co_await tt; // 1
    printf("main-3: %d\n", i);
    i = co_await tt; // 2
    printf("main-4: %d\n", i);
    i = co_await tt; // 3
    printf("main-5: %d\n", i);
    co_return 0;
}
```

```
main-1
test-1
main-2
main-3: 1
test-2
main-4: 2
test-3
main-5: 3
```

## Channels

Channels are modeled on golang; they are different from boost.asio channels in that they don't go through the executor.
Instead they directly context switch when possible.

```cpp
cobalt::promise<void> test(cobalt::channel<int> & chan)
{
  printf("Reader 1: %d\n", co_await chan.read());
  printf("Reader 2: %d\n", co_await chan.read());
  printf("Reader 3: %d\n", co_await chan.read());
}

cobalt::main co_main(int argc, char ** argv)
{
  cobalt::channel<int> chan{0u /* buffer size */};
  
  auto p = test(chan);
  
  printf("Writer 1\n");
  co_await chan.write(10);
  printf("Writer 2\n");
  co_await chan.write(11);
  printf("Writer 3\n");
  co_await chan.write(12);
  printf("Writer 4\n");
  
  co_await p;
  co_return 0u;
}
```

````
Writer-1
Reader-1: 10
Writer-2
Reader-1: 11
Writer-3
Reader-1: 12
Writer-4
````

## Ops

To make writing asio operations that have an early completion easier, cobalt has an op-helper:

```cpp
template<typename Timer>
struct wait_op : cobalt::op<system::error_code> // enable_op is to use ADL
{
  Timer & tim;

  wait_op(Timer & tim) : tim(tim) {}
  
  // this gets used to determine if it needs to suspend for the op
  void ready(cobalt::handler<system::error_code> h)
  {
    if (tim.expiry() < Timer::clock_type::now())
      h(system::error_code(asio::error::operation_aborted));
  }
  
  // this gets used to initiate the op if ti needs to suspend
  void initiate(cobalt::completion_handler<system::error_code> complete)
  {
    tim.async_wait(std::move(complete));
  }
};

cobalt::main co_main(int argc, char ** argv)
{
  cobalt::steady_timer tim{co_await cobalt::this_coro::executor}; // already expired
  co_await wait_op(tim); // will not suspend, since its ready
}

```


## race

[`race`](doc/reference/race.adoc) let's you await multiple awaitables at once. 

```cpp
cobalt::promise<void> delay(int ms)
{
    asio::steady_timer tim{co_await cobalt::this_coro::executor};
    tim.expires_after(std::chrono::milliseconds(ms));
    co_await tim.async_wait(cobalt::use_op);
}

cobalt::main co_main(int argc, char ** argv)
{
  auto res = co_await race(delay(100), delay(50));
  asert(res == 1); // delay(50) completes earlier, delay(100) is not cancelled  
  co_return 0u;
}
```

