# On the costs of asynchronous abstractions

The biggest force behind the evolution of
[Boost.Redis](https://github.com/boostorg/redis) was my struggling in
coming up with a high-level connection abstraction that was capable of
multiplexing Redis commands from independent sources while
concurrently handling server pushes. This journey taught me many
important lessons, many of which are related to the design and
performance of asynchronous programs based on Boost.Asio.

In this article I will share some of the lessons learned, specially
those related to the performance costs of _abstractions_ such as
`async_read_until` that tend to overschedule into the event-loop. In
this context I will also briefly comment on how the topics discussed
here influenced my views on the proposed
[P2300](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2023/p2300r7.html)
(a.k.a. Senders and Receivers), which is likely to become the basis of
networking in upcoming C++ standards.

Although the analysis presented in this article uses the Redis communication
protocol for illustration I expect it to be useful in general since
[RESP3](https://github.com/antirez/RESP3/blob/master/spec.md) shares
many similarities with other widely used protocols such as HTTP. 

## Parsing `\r\n`-delimited messages

The Redis server communicates with its clients by exchanging data
serialized in
[RESP3](https://github.com/antirez/RESP3/blob/master/spec.md) format.
Among the data types supported by this specification, the
`\r\n`-delimited messages are some of the most frequent in a typical
session. The table below shows some examples

  Command  | Response | Wire format   | RESP3 name
  ---------|----------|---------------|---------------------
  PING     |  PONG    | `+PONG\r\n`   | simple-string
  INCR     |  42      | `:42\r\n`     | number
  GET      |  null    | `_\r\n`       | null

Redis also supports command pipelines, which provide a way of
optimizing round-trip times by batching commands. A pipeline composed
by the commands shown in the previous table look like this

```
                       |   Sent in a     |
                       |   single write  |
+--------+             |                 |              +-------+
|        | -------->    PING + INCR + GET     --------> |       |
|        |                                              |       |
| Client |                                              | Redis |
|        |                                              |       |
|        | <-------- "+PONG\r\n:42\r\n_\r\n"  <-------- |       |
+--------+            |<------>|<---->|<-->|            +-------+
                      |                    |
                      |     Responses      |
```

Messages that use delimiters are so common in networking that a
facility called `async_read_until` for reading them incrementally from
a socket is already part of Boost.Asio. The coroutine below uses it to
print message contents to the screen

```cpp
awaitable<void> parse_resp3_simple_msgs(tcp::socket socket)
{
  for (std::string buffer;;) {
     auto n = co_await async_read_until(socket, dynamic_buffer(buffer), "\r\n");

     std::cout << buffer.substr(1, n - 3) << std::endl;

     // Consume the buffer.
     buffer.erase(0, n);
  }
}
```

If we pay attention to the buffer content as it is parsed by the code
above we can see it is rotated fairly often, for example

```
   "+PONG\r\n:100\r\n+OK\r\n_\r\n"
   ":100\r\n+OK\r\n_\r\n"
   "+OK\r\n_\r\n"
   "_\r\n"
   ""
```

When I first realized these, apparently excessive, buffer rotations I
was concerned they would impact the performance of Boost.Redis in a
severe way.  To measure the magnitude of this impact I came up with an
experimental implementation of Asio's `dynamic_buffer` that consumed
the buffer less eagerly than the `std::string::erase` function used
above. For that, the implementation increased a buffer offset up
to a certain threshold and only then triggered a (larger) rotation.
This is illustrated in the diagram below

```
    |<---- offset threshold ---->|
    |                            |
   "+PONG\r\n:100\r\n+OK\r\n_\r\n+PONG\r\n"
    |                                       # Initial offset

   "+PONG\r\n:100\r\n+OK\r\n_\r\n+PONG\r\n"
    |<------>|                              # After 1st message

   "+PONG\r\n:100\r\n+OK\r\n_\r\n+PONG\r\n"
    |<-------------->|                      # After 2nd message

   "+PONG\r\n:100\r\n+OK\r\n_\r\n+PONG\r\n"
    |<--------------------->|               # After 3rd message

   "+PONG\r\n:100\r\n+OK\r\n_\r\n+PONG\r\n"
    |<-------------------------->|          # Threshold crossed after the 4th message

   "+PONG\r\n"
    |                                       # After rotation
```

After comparing the performance differences between the two versions I
was surprised there wasn't any! But that was also very suspicious
since some RESP3 aggregate types contain a considerable number of
separators.  For example, a map with two pairs `[(key1, value1),
(key2, value2)]` encoded in RESP3 requires ten rotations in total

```
   "%2\r\n$4\r\nkey1\r\n$6\r\nvalue1\r\n$4\r\nkey2\r\n$6\r\nvalue2\r\n"
   "$4\r\nkey1\r\n$6\r\nvalue1\r\n$4\r\nkey2\r\n$6\r\nvalue2\r\n"
   "key1\r\n$6\r\nvalue1\r\n$4\r\nkey2\r\n$6\r\nvalue2\r\n"
   "$6\r\nvalue1\r\n$4\r\nkey2\r\n$6\r\nvalue2\r\n"
   ...
```

It was evident something more costly was shadowing the buffer
rotations. But it couldn't be the search for the separator since it
performs equivalently to rotations. It is also easy to show that the
overhead is not related to any IO operation since the problem persists
if the buffer is never consumed (which causes the function to be
called with the same string repeatedly).  Once these two factors
are removed from the table, we are driven into the conclusion that
calling `async_read_until` has an intrinsic cost, let us see what
that is.

### Async operations that complete synchronously considered harmful

Assume the scenario described earlier where `async_read_until` is used
to parse multiple `\r\n`-delimited messages. The following is a
detailed description of what happens behind the scenes

   1. `async_read_until` calls `socket.async_read_some` repeatedly
      until the separator `\r\n` shows up in the buffer

```
   "<read1>"                                # Read 1: needs more data.
   "<read1><read2>"                         # Read 2: needs more data.
   "<read1><read2>"                         # Read 3: needs more data.
   "<read1><read2><read3>"                  # Read 4: needs more data.
   "<read1><read2><read3>\r\n<bonus bytes>" # separator found, done.
```

   2. The last call to `socket.async_read_some` happens to read past
      the separator `\r\n` (depicted as `<bonus bytes>` above),
      resulting in bonus (maybe incomplete) messages in the buffer

```
    | 1st async_read_some  |         2nd async_read_some            |
    |                      |                                        |
    "+message content here \r\n:100\r\n+OK\r\n_\r\n+incomplete respo"
    |                          |                   |                |
    |      Message wanted      |<-- bonus msgs --->|<--incomplete-->|
    |                          |                         msg        |
    |                          |                                    |
    |                          |<---------- bonus bytes ----------->|
```

   3. The buffer is consumed and `async_read_until` is called again.
      However, since the buffer already contains the next message this
      is an IO-less call

```
   ":100\r\n+OK\r\n_\r\n+not enough byt"
    |                   |              |
    |  No IO required   | Need more    |
    |  to parse these   | data         |   
    |  messages.        |              |
```

The fact that step 3. doesn't perform any IO implies the operation can
complete synchronously, but because this is an asynchronous function
Boost.Asio by default won't call the continuation before the
function returns. The implementation must therefore enqueue it for
execution, as depicted below

```
   OP5 ---> OP4 ---> OP3 ---> OP2 ---> OP1    # Reschedules the continuation
                                        |
       OP1 schedules its continuation   |
    +-----------------------------------+
    |
    |
   OP6 ---> OP5 ---> OP4 ---> OP3 ---> OP2    # Reschedules the continuation
                                        |
       OP2 schedules its continuation   |
    +-----------------------------------+
    |
    |
   OP7 ---> OP6 ---> OP5 ---> OP4 ---> OP3    
```

When summed up, the excessive rescheduling of continuations lead to
performance degradation at scale. But since this is an event-loop
there is no way around rescheduling as doing otherwise would mean
allowing a task to monopolize the event-loop, preventing other tasks
from making progress. The best that can be done is to avoid
_overscheduling_, so let us determine how much rescheduling is too
much.

## The intrinsic latency of an event-loop 

An event-loop is a design pattern originally used to handle events
external to the application, such as GUIs, networking and other forms
of IO.  If we take this literally, it becomes evident that the way
`async_read_until` works is incompatible with an event-loop since
_searching for the separator_ is not an external event and as such
should not have to be enqueued for execution.

Once we constrain ourselves to events that have an external origin,
such as anything related to IO and including any form of IPC, the
scheduling overhead is reduced considerably since the latency
of the transport layer eclipses whatever time it takes to schedule the
continuation, for example, according to
[these](https://www.boost.org/doc/libs/develop/libs/cobalt/doc/html/index.html#posting_to_an_executor)
benchmarks, the time it takes to schedule a task in the
`asio::io_context ` is approximately `50ns`.

To give the reader an idea about the magnitude of this number, if
rescheduling alone were to account for 1% of the runtime of an app
that uses asynchronous IO to move around data in chunks of size 128kb,
then this app would have a throughput of approximately 24Gbs.  At such
high throughput multiple other factors kick in before any scheduling
overhead even starts to manifest.

It is therefore safe to say that only asynchronous operations that
don't perform or are not bound to any IO are ever likely to
overschedule in the sense described above. Those cases can be usually
avoided, this is what worked for Boost.Redis

  1. `async_read_until` was replaced with calls to
     `socket.async_read_some` and an incremental parser that does not
     do any IO.

  2. Channel `try_` functions are used to check if send and receive
     operations can be called without suspension. For example,
     `try_send` before `async_send` and `try_receive` before
     `async_receive` ([see also](https://github.com/chriskohlhoff/asio/commit/fe4fd7acf145335eeefdd19708483c46caeb45e5)
     `try_send_via_dispatch` for a more aggressive optimization).

  3. Coalescing of individual requests into a single payload to reduce
     the number of necessary writes on the socket, this is only
     possible because Redis supports pipelining (good protocols
     help!).

  4. Increased the socket read sizes to 4kb to reduce the number of
     reads (which is outweighed by the costs of rotating data in the
     buffer).

  5. Dropped the `resp3::async_read` abstraction. When I started
     developing Boost.Redis there was convincing precedent for having
     a `resp3::async_read` function to read complete RESP3 messages
     from a socket

     Name                                   | Description
     ---------------------------------------|-------------------
     `asio::ip::tcp::async_read`            | Reads `n` bytes from a stream.
     `beast::http::async_read`              | Reads a complete HTTP message.
     `beast::websocket::stream::async_read` | Reads a complete Websocket message.
     `redis::async_read`                    | Reads a complete RESP3 message.

     It turns out however that this function is also vulnerable to
     immediate completions since in command pipelines multiple
     responses show up in the buffer after a call to
     `socket.async_read_some`. When that happens each call to
     `resp3::async_read` is IO-less.
   
Sometimes it is not possible to avoid asynchronous operations that
complete synchronously, in the following sections we will see how to
favor throughput over fairness in Boost.Asio.

### Calling the continuation inline

In Boost.Asio it is possible to customize how an algorithm executes
the continuation when an immediate completion occurs, this includes
the ability of calling it inline, thereby avoiding the costs of
excessive rescheduling.  Here is how it works

```cpp
// (default) The continuation is enqueued for execution, regardless of
// whether it is immediate or not.
async_read_until(socket, buffer, "\r\n", continuation);

// Immediate completions are executed in exec2 (otherwise equal to the
// version above). The completion is called inline if exec2 is the
// same executor that is running the operation.
async_read_until(socket, buffer, "\r\n", bind_immediate_executor(exec2, completion));
```

To compare the performance of both cases I have written a small
function that calls `async_read_until` in a loop with a buffer that is
never consumed so that all completions are immediate.  The version
below uses the default behaviour

```cpp
void read_safe(tcp::socket& s, std::string& buffer)
{
   auto continuation = [&s, &buffer](auto ec, auto n)
   {
      read_safe(s, buffer); // Recursive call
   };

   // This won't cause stack exhaustion because the continuation is 
   // not called inline but posted in the event loop.
   async_read_until(s, dynamic_buffer(buffer), "\r\n", continuation);
}
```

To optimize away some of the rescheduling the version below uses the
`bind_immediate_executor` customization to call the continuation
reentrantly and then breaks the stack from time to time to avoid
exhausting it

```cpp
void read_reentrant(tcp::socket& s, std::string& buffer)
{
   auto cont = [&](auto, auto)
   {
      read_reentrant(s, buffer); // Recursive call
   };

   // Breaks the callstack after 16 inline calls.
   if (counter % 16 == 0) {
      post(s.get_executor(), [cont](){cont({}, 0);});
      return;
   }

   // Continuation called reentrantly.
   async_read_until(s, dynamic_buffer(buffer), "\r\n",
      bind_immediate_executor(s.get_executor(), cont));
}
```

The diagram below shows what the reentrant chain of calls in the code
above look like from the event-loop point of view

```
   OP5 ---> OP4 ---> OP3 ---> OP2 ---> OP1a    # Completes immediately
                                        |
                                        |      
                 ...                    |
                                       OP1b    # Completes immediately
                                        |
           Waiting for OP5 to           |      
             reschedule its             |
             continuation              OP1c    # Completes immediately
                                        |
                                        |      
                 ...                    |
                                       OP1d    # Break the call-stack
                                        |
    +-----------------------------------+     
    |
   OP6 ---> OP5 ---> OP4 ---> OP3 ---> OP2
```

Unsurprisingly, the reentrant code is 3x faster than the one that
relies on the default behaviour (don't forget that this is a best case
scenario, in the general case not all completions are immediate).
Although faster, this strategy has some downsides

   - The overall operation is not as fast as possible since it still
     has to reschedule from time to time to break the call stack.  The
     less it reschedules the higher the risk of exhausting it.

   - It is too easy to forget to break the stack. For example, the
     programmer might decide to branch somewhere into another chain of
     asynchronous calls that also use this strategy.  To avoid
     exhaustion all such branches would have to be safeguarded with a
     manual rescheduling i.e. `post`.

   - Requires additional layers of complexity such as
     `bind_immediate_executor` in addition to `bind_executor`.
   
  - Non-compliat with more strict
    [guidelines](https://en.wikipedia.org/wiki/The_Power_of_10:_Rules_for_Developing_Safety-Critical_Code)
    that prohibits reentrat code.

   - There is no simple way of choosing the maximum allowed number of
     reentrant calls for each function in a way that covers different
     use cases and users. Library writers and users would be tempted
     into using a small value reducing the performance advantage.

   - If the socket is always ready for reading the task will
     monopolize IO for up to `16` interactions which might cause
     stutter in unrelated tasks as depicted below

```
                                 Unfairness

          +----+----+----+    +----+----+----+    +----+----+----+
Socket-1  |    |    |    |    |    |    |    |    |    |    |    |
          +----+----+----+----+----+----+----+----+----+----+----+----+     
Socket-2                 |    |              |    |              |    | 
                         +----+              +----+              +----+
```

From the aesthetic point of view the code above is also unpleasant as
it breaks the function asynchronous contract by injecting a reentrant
behaviour. It gives me the same kind of feeling I have about
[recursive
mutexes](http://www.zaval.org/resources/library/butenhof1.html).

Note: It is worth mentioning here that a similar
[strategy](https://github.com/NVIDIA/stdexec/blob/6f23dd5b1d523541ce28af32fc2603403ebd36ed/include/exec/trampoline_scheduler.hpp#L52)
is used to break the call stack of repeating algorithms in
[stdexec](https://github.com/NVIDIA/stdexec), but in this time
based on
[P2300](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2023/p2300r7.html)
and not on Boost.Asio.

### Coroutine tail-calls

In the previous section we have seen how to avoid overscheduling by
instructing the asynchronous operation to call the completion inline
on immediate completion.  It turns out however that coroutine support
for _tail-calls_ provide a way to completely sidestep this problem.
This feature is described by
[Lewis Baker](https://lewissbaker.github.io/2020/05/11/understanding_symmetric_transfer)
as follows

> A tail-call is one where the current stack-frame is popped before
> the call and the current function’s return address becomes the
> return-address for the callee. ie. the callee will return directly
> the the [sic] caller of this function.

This means (at least in principle) that a library capable of using
tail-calls when an immediate completion occurs neither has to
reschedule the continuation nor call it inline. To test how this
feature compares to the other styles I have used Boost.Cobalt. The
code looks as follows

```cpp
// Warning: risks unfairness and starvation of other tasks.
task<void> read_until_unfair()
{
    for (int i = 0; i != repeat; ++i) {
        co_await async_read_until(s, dynamic_buffer(buffer), "\r\n", cobalt::use_op);
    }
}
```

The result of this comparison as listed in the table below

Time/s | Style     | Configuration               | Library     
-------|-----------|-----------------------------|-------------
 1,0   | Coroutine | `await_ready` optimization  | Boost.Cobalt
 4.8   | Callback  | Reentant                    | Boost.Asio  
10.3   | Coroutine | `use_op`                    | Boost.Cobalt
14.9   | Callback  | Regular                     | Boost.Asio  
15.6   | Coroutine | `asio::deferred`            | Boost.Asio  

As the reader can see, `cobalt::use_op` ranks 3rd and is considerably
faster (10.3 vs 15.6) than the Asio equivalent that uses
default-rescheduling. However, by trading rescheduling with tail-calls
the code above can now monopolize the event-loop, resulting in
unfairness if the socket happens to receive data at a higher rate
than other tasks. If by chance data is received continuously
on a socket that is always ready for reading, other tasks will starve

```
                                 Starvation

          +----+----+----+----+----+----+----+----+----+----+----+----+
Socket-1  |    |    |    |    |    |    |    |    |    |    |    |    |
          +----+----+----+----+----+----+----+----+----+----+----+----+     

Socket-2          Starving ...                                          
                                                                       
```

To avoid this problem the programmer is forced to reschedule from time
to time, in the same way we did for the reentrant calls

```cpp
task<void> read_until_fair()
{
    for (int i = 0; i != repeat; ++i) {
        if (repeat % 16 == 0) {
            // Reschedules to address unfairness and starvation of 
            // other tasks.
            co_await post(cobalt::use_op);
            continue;
        }
        
        co_await async_read_until(s, dynamic_buffer(buffer), "\r\n", cobalt::use_op);
    }
}
```

Delegating fairness-safety to applications is a dangerous game.
This is a
[problem](https://tokio.rs/blog/2020-04-preemption) the Tokio
community had to deal with before Tokio runtime started enforcing
rescheduling (after 256 successful operations)

> If data is received faster than it can be processed, it is possible
> that more data will have already been received by the time the
> processing of a data chunk completes. In this case, .await will
> never yield control back to the scheduler, other tasks will not be
> scheduled, resulting in starvation and large latency variance.

> Currently, the answer to this problem is that the user of Tokio is
> responsible for adding yield points in both the application and
> libraries. In practice, very few actually do this and end up being
> vulnerable to this sort of problem.

### Safety in P2300 (Senders and Receivers)

As of this writing, the C++ standards committee (WG21) has been
pursuing the standardization of a networking library for almost 20
years. One of the biggest obstacles that prevented it from happening
was a disagreement on what the _asynchronous model_ that underlies
networking should look like. Until 2021 that model was basically
Boost.Asio _executors_,  but in this
[poll](https://www.reddit.com/r/cpp/comments/q6tgod/c_committee_polling_results_for_asynchronous/)
the committee decided to abandon that front and concentrate efforts on
the new [P2300](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2023/p2300r7.html)
proposal, also known as _senders and receivers_. The decision was
quite [abrupt](https://isocpp.org/files/papers/P2464R0.html)

> The original plan about a week earlier than the actual writing of
> this paper was to write a paper that makes a case for standardizing
> the Networking TS.

and opinions turned out to be very strong against Boost.Asio (see
[this](https://api.csswg.org/bikeshed/?force=1&url=https://raw.githubusercontent.com/brycelelbach/wg21_p2459_2022_january_library_evolution_poll_outcomes/main/2022_january_library_evolution_poll_outcomes.bs)
for how each voter backed their vote)

> The whole concept is completely useless, there's no composed code
> you can write with it.

The part of that debate that interests us most here is stated in
[P2471](https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2021/p2471r1.pdf),
that compares Boost.Asio with P2300

> Yes, default rescheduling each operation and default not
> rescheduling each operation, is a poor trade off. IMO both options
> are poor. The one good option that I know of that can prevent stack
> exhaustion is first-class tail-recursion in library or language

> ASIO has chosen to require that every async operation must schedule
> the completion on a scheduler (every read, every write, etc..).

> sender/receiver has not decided to
> require that the completion be scheduled. 

> This is why I consider tail-call the only good solution. Scheduling
> solutions are all inferior (give thanks to Lewis for this shift in
> my understanding :) ).

Although tail-calls solve the problem of stack-exhaustion as we have
seen above, it makes the code vulnerable to unfairness and starvation
and therefore it is not an alternative to default-rescheduling as the
quotation above is implying.  To deal with the lack of
default-rescheduling, libraries and applications built on top of P2300
have to address the aforementioned problems, layer after layer. For
example,
[stdexec](https://github.com/NVIDIA/stdexec) has invented something
called
_[trampoline-scheduler](https://github.com/NVIDIA/stdexec/blob/e7cd275273525dbc693f4bf5f6dc4d4181b639e4/include/exec/trampoline_scheduler.hpp)_
to protect repeating algorithms such as `repeat_effect_until` from
exhausting the stack. This construct however is built around
reentracy, allowing
[sixteen](https://github.com/NVIDIA/stdexec/blob/83cdb92d316e8b3bca1357e2cf49fc39e9bed403/include/exec/trampoline_scheduler.hpp#L52)
levels of inline calls by default. While in Boost.Asio it is possible to use
reentracy as an optimization for a corner cases, here it is made its
_modus operandi_, the downsides of this approach have already been stated in a
previous section so I won't repeat it here.

Also the fact that a special scheduler is needed by specific
algorithms is a problem on its own since it contradicts one of the
main selling points of P2300 which is that of being _generic_. For
example, [P2464R0](https://isocpp.org/files/papers/P2464R0.html) uses
the code below as an example

```cpp
void
run_that_io_operation(
    scheduler auto sched,
    sender_of<network_buffer> auto wrapping_continuation)
{
   // snip
}
```

and states

> I have no idea what the sched's concrete type is. I have no idea
> what the wrapping_continuation's concrete type is. They're none of
> my business, ...

Hence, by being generic, the algorithms built on top of P2300 are also
unsafe (against stack-exhaustion, unfairness and starvation).  Otherwise,
if library writers require a specific scheduler to ensure safety, then
the algorithms become automatically non-generic, pick your poison!

The proposers of P2300 claim that it doesn't address safety because it
should be seen as the low-level building blocks of asynchronous
programming and that its the role of higher-level libraries, to deal
with that. This claim however does not hold since, as we have just
seen, Boost.Asio also provides those building blocks but does so in a
safe way.  In fact during the whole development of Boost.Redis I never
had to think about these kinds of problems because safety is built
from the ground up.

### Avoiding coroutine suspension with `await_ready`

Now let us get back to the first place in the table above, which uses
the `await_ready` optimization from Boost.Cobalt. This API provides
users with the ability to avoid coroutine suspension altogether in
case the separator is already present in the buffer.  It works by
defining a `struct` with the following interface

```cpp
struct read_until : cobalt::op<error_code, std::size_t> {
   ...

   void ready(cobalt::handler<error_code, std::size_t> handler) override
   {
      // Search for the separator in buffer and call the handler if found
   }

   void initiate(cobalt::completion_handler<error_code, std::size_t> complete) override
   {
      // Regular call to async_read_until.
      async_read_until(socket, buffer, delim, std::move(complete));
   }
};
```

and the code that uses it

```cpp
for (int i = 0; i != repeat; ++i) {
    co_await read_until(socket, dynamic_buffer(buffer));
}
```

In essence, what the code above does is to skip a call to
`async_read_unil` by first checking with the ready function whether
the forthcoming operation is going to complete immediately.  The
nice thing about it is that the programmer can use this optimization
only when a performance bottleneck is detected, without planing for it
in advance.  The drawback however is that it requires reimplementing
the search for the separator in the body of the `ready` function,
defeating the purpose of using `async_read_until` in first place as
(again) it would have been simpler to reformulate the operation in
terms of `socket.async_read_some` directly.

## Acknowledgements

Thanks to Klemens Morgenstern for answering questions about
Boost.Cobalt.

