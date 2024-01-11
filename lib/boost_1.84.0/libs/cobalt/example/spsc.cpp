//
// Copyright (c) 2023 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt.hpp>

#include <boost/asio/post.hpp>
#include <boost/lockfree/spsc_queue.hpp>

namespace cobalt = boost::cobalt;

/// This is a simple class making a lockfree::spsc_queue awaitable. It's not movable, since the spsc_queue isn't.
template<typename T, typename ... Options>
struct awaitable_spsc_queue
{
  boost::lockfree::spsc_queue<T, Options...> queue;

  template<typename ... Args>
  awaitable_spsc_queue(Args && ... args) : queue(std::forward<Args>(args)...) {}

  // if the queue gets destroyed, destroying the awaiters is all we can do.
  ~awaitable_spsc_queue()
  {
    if (auto r = reader.load(); r != nullptr)
      std::coroutine_handle<void>::from_address(r).destroy();

    if (auto w = writer.load(); w != nullptr)
      std::coroutine_handle<void>::from_address(w).destroy();
  }
  // to avoid locks, put the coroutine handles into atomics.
  std::atomic<void*> reader{nullptr}, writer{nullptr};

  // capture the read & write executor
  cobalt::executor read_executor, write_executor;

  // the awaitable to read a value from the queue
  struct read_op
  {
    awaitable_spsc_queue * this_;

    // if reads are available we don't need to suspend
    bool await_ready() {return this_->queue.read_available();}

    // The suspend implementation. We expect the coroutine promise to have an associated executor.
    // If it doesn't resumption will end up on the system_executor
    template<typename Promise>
    void await_suspend(std::coroutine_handle<Promise> h)
    {
      // Capture the read_executor.
      if constexpr (requires {h.promise().get_executor();})
        this_->read_executor = h.promise().get_executor();
      else
        this_->read_executor = cobalt::this_thread::get_executor();

      // Make sure there's only one coroutine awaiting the read
      assert(this_->reader == nullptr);
      // Store the handle of the awaiter.
      this_->reader.store(h.address());
    }
    T await_resume()
    {
      T res;
      // Grab the value from the queue
      this_->queue.pop(res);
      // if a writer is waiting post it to complete on it's thread
      auto w = cobalt::unique_handle<void>::from_address(this_->writer.exchange(nullptr));
      if (w)
        boost::asio::post(this_->write_executor, std::move(w));

      return res;
    }
  };

  // The function to be used with for reading.
  read_op read() {return {this};}

  // the awaitable to write a value from the queue
  struct write_op
  {
    awaitable_spsc_queue * this_;
    // The value to write
    T value;

    // If there is room in the queue we don't need to suspend
    bool await_ready() {return this_->queue.write_available();}

    // The suspend implementation. Similar to the writer one.
    template<typename Promise>
    void await_suspend(std::coroutine_handle<Promise> h)
    {
      // Capture the read_executor.
      if constexpr (requires {h.promise().get_executor();})
        this_->write_executor = h.promise().get_executor();
      else
        this_->write_executor = cobalt::this_thread::get_executor();

      assert(this_->writer == nullptr);
      this_->writer.store(h.address());
    }

    // Resume the coroutine, this is where we do the actual write.
    // The reason is that we suspend when we cannot write because the queue is full.
    void await_resume()
    {
      this_->queue.push(std::move(value));
      // if a writer is waiting post it
      auto r = cobalt::unique_handle<void>::from_address(this_->reader.exchange(nullptr));
      if (r)
        boost::asio::post(this_->read_executor, std::move(r));
    }
  };

  // The function to be used for writing.
  write_op write(T value) {return {this, std::move(value)};}
};

// Dummy thread blasting out values.
cobalt::thread thr(awaitable_spsc_queue<int> & q)
{
  for (int i = 0; i <= 100000000; i++)
    co_await q.write(i);

  co_await q.write(-1);
}

cobalt::main co_main(int argc, char * argv[])
{
  awaitable_spsc_queue<int> queue{1024};
  auto t = thr(queue);

  // dummy consumer
  for (auto val = co_await queue.read();
       val >= 0;
       val = co_await queue.read());

  co_await t;
  co_return 0;
}