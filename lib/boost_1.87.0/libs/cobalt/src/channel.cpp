//
// Copyright (c) 2022 Klemens Morgenstern (klemens.morgenstern@gmx.net)
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/cobalt/channel.hpp>
#include <boost/asio/defer.hpp>

namespace boost::cobalt
{

channel<void>::~channel()
{
  while (!read_queue_.empty())
    read_queue_.front().awaited_from.reset();

  while (!write_queue_.empty())
    write_queue_.front().awaited_from.reset();

}

void channel<void>::close()
{
  is_closed_ = true;
  while (!read_queue_.empty())
  {
    auto & op = read_queue_.front();
    op.unlink();
    op.cancelled = true;
    op.cancel_slot.clear();
    if (op.awaited_from)
      asio::defer(executor_, std::move(op.awaited_from));
  }
  while (!write_queue_.empty())
  {
    auto & op = write_queue_.front();
    op.unlink();
    op.cancelled = true;
    op.cancel_slot.clear();
    if (op.awaited_from)
      asio::defer(executor_, std::move(op.awaited_from));
  }
}

system::result<void>  channel<void>::read_op::await_resume(const struct as_result_tag &)
{
  if (cancel_slot.is_connected())
    cancel_slot.clear();

  if (cancelled)
    return {system::in_place_error, asio::error::operation_aborted};

  if (!direct)
    chn->n_--;
  if (!chn->write_queue_.empty())
  {
    auto &op = chn->write_queue_.front();
    BOOST_ASSERT(chn->read_queue_.empty());
    if (op.await_ready())
    {
      op.unlink();
      BOOST_ASSERT(op.awaited_from);
      asio::post(
          chn->executor_, std::move(op.awaited_from));
    }
  }
  return {system::in_place_value};
}

void channel<void>::read_op::await_resume()
{
  await_resume(as_result_tag{}).value(loc);
}

std::tuple<system::error_code> channel<void>::read_op::await_resume(const struct as_tuple_tag & )
{
  return await_resume(as_result_tag{}).error();
}


system::result<void> channel<void>::write_op::await_resume(const struct as_result_tag &)
{
  if (cancel_slot.is_connected())
    cancel_slot.clear();
  if (cancelled)
    return {system::in_place_error, asio::error::operation_aborted};
  if (!direct)
    chn->n_++;

  if (!chn->read_queue_.empty())
  {
    auto & op = chn->read_queue_.front();
    BOOST_ASSERT(chn->write_queue_.empty());
    if (op.await_ready())
    {
      op.unlink();
      BOOST_ASSERT(op.awaited_from);
      asio::post(
          chn->executor_, std::move(op.awaited_from));
    }
  }
  return {system::in_place_value};
}


void channel<void>::write_op::await_resume()
{
  await_resume(as_result_tag{}).value(loc);
}


std::tuple<system::error_code> channel<void>::write_op::await_resume(const struct as_tuple_tag & )
{
  return await_resume(as_result_tag{}).error();
}

}
