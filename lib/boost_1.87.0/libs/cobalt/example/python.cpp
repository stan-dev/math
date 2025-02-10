// Copyright (c) 2023 Klemens D. Morgenstern
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <nanobind/nanobind.h>
#include <boost/asio/executor_work_guard.hpp>
#include <boost/asio/steady_timer.hpp>

#include <boost/cobalt.hpp>
#include <thread>

/** In this example we'll connect cobalt and
 * pythons asyncio using nanobind (a C++17 successor to pybind11).
 *
 */

using namespace boost;
namespace py = nanobind;


// Small helper function to get the current event loop for python
py::object get_loop()
{
  auto mod = py::module_::import_("asyncio");
  auto getter = mod.attr("get_event_loop");
  return getter();
}

// An asio::any_io_executor compatible wrapper around the event loop.
// all the query functions are for that compatibility.
struct python_executor
{
  python_executor(py::object loop = ::get_loop()) noexcept
      : m_ptr(loop.ptr())
  {
    // the asio::execution_context needs to present, we put it into the event loop.
    if (!py::hasattr(loop, "__asio_execution_context"))
    {
      auto ptr = std::make_unique<asio::execution_context>();
      context_ = ptr.get();
      py::setattr(loop, "__asio_execution_context",
                  py::cast(std::move(ptr)));
    }
    else
      context_ = py::cast<std::unique_ptr<asio::execution_context>&>(
                          py::getattr(loop, "__asio_execution_context")
                      ).get();
  }

  asio::execution_context &query(asio::execution::context_t) const 
  {
    return *context_;
  }

  static constexpr asio::execution::blocking_t
  query(asio::execution::blocking_t) noexcept
  {
    return asio::execution::blocking.never;
  }

  // this function takes the function F and runs it on the event loop.
  template<class F>
  void
  execute(F f) const
  {
    py::gil_scoped_acquire lock;
    py::handle loop(m_ptr);
    struct wrapper // override the const.
    {
      mutable F impl;

      void operator()( ) const
      {
        std::move(impl)();
      }
    };
    loop.attr("call_soon_threadsafe")(py::cpp_function(wrapper{std::move(f)}));
  }

  bool
  operator==(python_executor const &other) const noexcept
  {
    return m_ptr == other.m_ptr;
  }

  bool
  operator!=(python_executor const &other) const noexcept
  {
    return !(*this == other);
  }

  private:
  PyObject *m_ptr;
  asio::execution_context *context_;
};

// helper function so we can capture the currently active C++ exception into a python object.
py::object translate_current_exception()
{
  py::object locals = py::dict();
  locals["rethrow_"] = py::cpp_function([]{throw;});
  PyRun_String(R"(
try:
    rethrow_();
    result = None
except Exception as e:
    result = e;

)", Py_single_input, PyEval_GetGlobals(), locals.ptr());
  return locals["result"];
}

struct py_coroutine
{
  struct promise_type
  {
    constexpr static std::suspend_always initial_suspend() noexcept { return {};}
    std::suspend_never final_suspend() noexcept
    {
      if (!done_)
        return_value(py::none());
      return {};
    }

    void unhandled_exception()
    {
      done_ = true;
      loop.attr("call_soon")(
          future.attr("set_exception"), translate_current_exception());
    }

    // cast errors will be caught by unhandled_exception.
    template<typename T>
    void return_value(T && value)
    {
      done_ = true;
      loop.attr("call_soon")(
          future.attr("set_result"), py::cast(std::forward<T>(value)));
    }

    template<typename T>
    std::suspend_always yield_value(T && value)
    {
      loop.attr("call_soon")(
          future.attr("set_result"), py::cast(std::forward<T>(value)));
      owner->handle_.reset(this);
      return {};
    }

    py::object loop;
    py::object future;

    py_coroutine get_return_object()
    {
      return py_coroutine{this};
    }

    bool done_ = false;

    using executor_type = asio::any_io_executor;
    executor_type exec_;
    const executor_type & get_executor() { return exec_;}

    py_coroutine * owner;
  };

  void initiate(py::object loop, py::object future)
  {
    if (!handle_)
      throw std::invalid_argument("Awaited invalid coroutine");

    if (handle_->done_)
      throw py::stop_iteration("coroutine completed");

    handle_->loop = std::move(loop);
    handle_->future = std::move(future);
    handle_->exec_ = python_executor(handle_->loop);
    handle_->owner = this;

    if (!cobalt::this_thread::has_executor())
      cobalt::this_thread::set_executor(handle_->exec_);
    std::coroutine_handle<promise_type>::from_promise(*handle_.release()).resume();
  }

  bool done() const
  {
    return !handle_ || handle_->done_;
  }

  py_coroutine(const py_coroutine & ) = delete;
  py_coroutine(py_coroutine && ) = default;

 private:
  explicit py_coroutine(promise_type * pro) : handle_(pro)
  {
  }
  struct deleter
  {
    void operator()(promise_type * p)
    {
      std::coroutine_handle<promise_type>::from_promise(*p).destroy();
    }
  };
  std::unique_ptr<promise_type, deleter> handle_;
};


struct py_awaitable
{
  py::object target;
  py_awaitable(py::object target) : target(std::move(target)) {}

  py::object result;
  py::object exception;
  constexpr bool await_ready() noexcept {return false;}

  template<typename T>
  bool await_suspend(std::coroutine_handle<T> h)
  {
    asio::any_io_executor exec;
    if constexpr (requires (T & p) {p.get_executor();})
      exec = h.promise().get_executor();
    else
      exec = cobalt::this_thread::get_executor();

    auto loop = get_loop();
    auto task = getattr(loop, "create_task")(target);

    struct wrapper
    {
      asio::any_io_executor exec;
      mutable cobalt::unique_handle<void> awaiter;
      py_awaitable * res;

      void operator()(py::object t) const
      {
        res->extract_result(t);
        asio::dispatch(exec, std::move(awaiter));
      }
    };

    if (py::cast<bool>(getattr(task, "done")()))
    {
      extract_result(task);
      return false;
    }

    getattr(task, "add_done_callback")(py::cpp_function(wrapper{
          std::move(exec),
          cobalt::unique_handle<void>{h.address()},
          this
        }));

    if constexpr (requires (T & p) {p.get_cancellation_slot();})
      h.promise().get_cancellation_slot().
      assign([c = getattr(task, "cancel")](asio::cancellation_type ct) {
        py::gil_scoped_acquire lock;
        c();
      });

    return true;
  }

  void extract_result(py::object task)
  {
    exception = getattr(task, "exception")();
    if (exception.is_none())
      result = getattr(task, "result")();
  }

  py::object await_resume()
  {
    if (!exception.is_none())
    {
      py::object locals = py::dict();
      PyRun_String(
          R"(def __rethrow__(ex):
    raise ex;
)",
          Py_single_input, PyEval_GetGlobals(), locals.ptr());
      locals["__rethrow__"](exception);
    }
    return std::move(result);
  }
};

cobalt::generator<int> test_generator()
{
  for (auto i = 1; i < 10; i++)
    co_yield i;
  co_return 10;
}

cobalt::promise<int> test_promise()
{
  asio::steady_timer tim{co_await cobalt::this_coro::executor,
                         std::chrono::milliseconds(100)};

  co_await tim.async_wait(cobalt::use_op);
  co_return 42;
}


cobalt::promise<void> await_py_coroutine(py_awaitable aw)
{
  auto res = co_await std::move(aw);
  printf("Python coroutine gave %s\n", py::str(res).c_str());
}


NB_MODULE(boost_cobalt_example_python, m)
{
  namespace execution = asio::execution;
  m.def("__rethrow_exception", &std::rethrow_exception);

  py::class_<std::unique_ptr<asio::execution_context>>(m, "__asio__execution_context");

  // use some inlined python to get a future
  py::object locals = py::dict();
  // language=python
  PyRun_String(
      R"(def __await_impl(self):
    import asyncio
    lp = asyncio.get_event_loop()
    ft = lp.create_future()
    self.initiate(lp, ft)
    res = yield from ft
    return res)", Py_single_input, PyEval_GetGlobals(), locals.ptr());
  py::object await_impl = locals["__await_impl"];

  // language=python
  PyRun_String(
      R"(async def __aiter_impl(self):
    while not self.done:
        yield await self
)", Py_single_input, PyEval_GetGlobals(), locals.ptr());

  py::object aiter_impl = locals["__aiter_impl"];

  py::class_<py_coroutine> ct(m, "__cobalt_coroutine");
  ct.def("initiate", &py_coroutine::initiate)
    .def_prop_ro("done", &py_coroutine::done);
  setattr(ct, "__await__", await_impl);
  setattr(ct, "__aiter__", aiter_impl);

  m.def("test_generator",
        []() -> py_coroutine
        {
            BOOST_COBALT_FOR(auto v, test_generator())
              co_yield v;
            co_return py::none();
        });

  m.def("test_promise",
        []() -> py_coroutine
        {
          co_return co_await test_promise();
        });
  m.def("test_py_promise",
        [](py::object obj) -> py_coroutine
        {
          co_await await_py_coroutine(std::move(obj));
          co_return py::none();
        });
}

