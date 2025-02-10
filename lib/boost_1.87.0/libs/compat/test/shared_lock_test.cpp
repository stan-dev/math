// Copyright 2023 Peter Dimov.
// Copyright 2023 Christian Mazakas.
// Copyright 2023 Ed Catmur
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

#ifdef _MSC_VER
#pragma warning(disable: 4530) // C++ exception handler used, but unwing semantics not enabled
#pragma warning(disable: 4577) // noexcept used with no exception handling mode specified
#endif

#include <boost/compat/shared_lock.hpp>

#include <boost/core/lightweight_test.hpp>
#include <boost/core/lightweight_test_trait.hpp>
#include <boost/throw_exception.hpp>
#include <boost/config.hpp>

#include <atomic>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <stdexcept>
#include <cstdio>

#define STATIC_ASSERT( ... ) static_assert( __VA_ARGS__, #__VA_ARGS__ )

#ifdef __clang__
#pragma clang diagnostic ignored "-Wself-move"
#endif
#if defined(__GNUC__) && __GNUC__ >= 13
#pragma GCC diagnostic ignored "-Wself-move"
#endif

struct invalid_lock_use: public std::runtime_error
{
    invalid_lock_use(): std::runtime_error( "Invalid lock use" )
    {
    }
};

struct dummy_lock {
private:
  int lock_shared_count_ = 0;
  int unlock_shared_count_ = 0;

  int lock_unique_count_ = 0;
  int unlock_unique_count_ = 0;

public:
  dummy_lock() = default;
  ~dummy_lock() {
    BOOST_TEST_EQ( lock_shared_count_, unlock_shared_count_ );
    BOOST_TEST_EQ( lock_unique_count_, unlock_unique_count_ );
  }

  void lock() {
    if ( lock_shared_count_ != unlock_shared_count_ ) {
      boost::throw_exception( invalid_lock_use(), BOOST_CURRENT_LOCATION );
    }
    ++lock_unique_count_;
  }

  void unlock() { ++unlock_unique_count_; }

  bool try_lock() {
    if ( lock_shared_count_ != unlock_shared_count_ ) {
      return false;
    }

    ++lock_unique_count_;
    return true;
  }

  void lock_shared() {
    if ( lock_unique_count_ != unlock_unique_count_ ) {
      boost::throw_exception( invalid_lock_use(), BOOST_CURRENT_LOCATION );
    }
    ++lock_shared_count_;
  }

  bool try_lock_shared() {
    if ( lock_unique_count_ != unlock_unique_count_ ) {
      return false;
    }
    ++lock_shared_count_;
    return true;
  }

  void unlock_shared() { ++unlock_shared_count_; }

  int shared_lock_count() const noexcept { return lock_shared_count_; }
};

namespace {

using shared_lock_type = boost::compat::shared_lock<dummy_lock>;

// verify that our dummy_lock throws when we hope it does
void sanity_tests() {
  dummy_lock sp;
  sp.lock();
  BOOST_TEST_THROWS( sp.lock_shared(), invalid_lock_use );
  sp.unlock();

  sp.lock_shared();
  BOOST_TEST_THROWS( sp.lock(), invalid_lock_use );
  sp.unlock_shared();
}

void default_constructor() {
  // shared_lock() noexcept;

  BOOST_TEST_TRAIT_SAME( typename shared_lock_type::mutex_type, dummy_lock );
  BOOST_TEST_TRAIT_TRUE(
      (std::is_nothrow_default_constructible<shared_lock_type>));

  shared_lock_type lock;
  BOOST_TEST_EQ( lock.mutex(), nullptr );
  BOOST_TEST( !lock.owns_lock() );
  BOOST_TEST( !lock );
}

void locking_construtor() {
  // explicit shared_lock( mutex_type& m );

  dummy_lock sp;

  {
    shared_lock_type lock( sp );
    BOOST_TEST_EQ( lock.mutex(), std::addressof( sp ) );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );
    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  }
}

void deferred_constructor() {
  // shared_lock( mutex_type& m, defer_lock_t t ) noexcept;

  dummy_lock sp;

  shared_lock_type lock( sp, std::defer_lock );

  BOOST_TEST_EQ( lock.mutex(), std::addressof( sp ) );
  BOOST_TEST( !lock.owns_lock() );
  BOOST_TEST( !lock );
  BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
}

void try_lock_constructor() {
  // shared_lock( mutex_type& m, try_to_lock_t t );

  dummy_lock sp;

  {
    shared_lock_type lock( sp, std::try_to_lock );
    BOOST_TEST_EQ( lock.mutex(), std::addressof( sp ) );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );
    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  }

  {
    sp.lock();

    shared_lock_type lock( sp, std::try_to_lock );
    BOOST_TEST_EQ( lock.mutex(), std::addressof( sp ) );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );
    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );

    sp.unlock();
  }
}

void adopt_lock_constructor() {
  // shared_lock(mutex_type& m, adopt_lock_t);

  dummy_lock sp;
  sp.lock_shared();

  {
    shared_lock_type lock( sp, std::adopt_lock );
    BOOST_TEST_EQ( lock.mutex(), std::addressof( sp ) );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );
    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  }

  sp.lock();
  sp.unlock();
}

void move_constructor() {
  BOOST_TEST_TRAIT_TRUE(
      (std::is_nothrow_move_constructible<shared_lock_type>));

  {
    dummy_lock sp;

    shared_lock_type lock( sp );
    shared_lock_type lock2( std::move( lock ) );
    BOOST_TEST_EQ( lock2.mutex(), std::addressof( sp ) );
    BOOST_TEST( lock2.owns_lock() );
    BOOST_TEST( lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  }

  {
    dummy_lock sp;

    shared_lock_type lock( sp, std::defer_lock );
    shared_lock_type lock2( std::move( lock ) );
    BOOST_TEST_EQ( lock2.mutex(), std::addressof( sp ) );
    BOOST_TEST( !lock2.owns_lock() );
    BOOST_TEST( !lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
  }

  {
    dummy_lock sp;

    shared_lock_type lock;
    shared_lock_type lock2( std::move( lock ) );
    BOOST_TEST_EQ( lock2.mutex(), nullptr );
    BOOST_TEST( !lock2.owns_lock() );
    BOOST_TEST( !lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
  }
}

void move_assignment() {
  BOOST_TEST_TRAIT_TRUE( (std::is_nothrow_move_assignable<shared_lock_type>));

  {
    // lhs not bound, not locked
    // rhs not bound, not locked
    dummy_lock sp, sp2;

    shared_lock_type lock;
    shared_lock_type lock2;

    lock2 = std::move( lock );
    BOOST_TEST_EQ( lock2.mutex(), nullptr );
    BOOST_TEST( !lock2.owns_lock() );
    BOOST_TEST( !lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
    BOOST_TEST_EQ( sp2.shared_lock_count(), 0 );
  }

  {
    // lhs not bound, not locked
    // rhs bound, locked
    dummy_lock sp;

    shared_lock_type lock( sp );
    shared_lock_type lock2;

    lock2 = std::move( lock );
    BOOST_TEST_EQ( lock2.mutex(), std::addressof( sp ) );
    BOOST_TEST( lock2.owns_lock() );
    BOOST_TEST( lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  }

  {
    // lhs not bound, not locked
    // rhs bound, not locked
    dummy_lock sp, sp2;

    shared_lock_type lock( sp, std::defer_lock );
    shared_lock_type lock2;

    lock2 = std::move( lock );
    BOOST_TEST_EQ( lock2.mutex(), std::addressof( sp ) );
    BOOST_TEST( !lock2.owns_lock() );
    BOOST_TEST( !lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
    BOOST_TEST_EQ( sp2.shared_lock_count(), 0 );
  }

  {
    // lhs bound, locked
    // rhs not bound, not locked

    dummy_lock sp, sp2;

    shared_lock_type lock;
    shared_lock_type lock2( sp2 );

    lock2 = std::move( lock );

    BOOST_TEST_EQ( lock2.mutex(), nullptr );
    BOOST_TEST( !lock2.owns_lock() );
    BOOST_TEST( !lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
    BOOST_TEST_EQ( sp2.shared_lock_count(), 1 );
  }

  {
    // lhs bound, not locked
    // rhs not bound, not locked

    dummy_lock sp, sp2;

    shared_lock_type lock;
    shared_lock_type lock2( sp, std::defer_lock );

    lock2 = std::move( lock );

    BOOST_TEST_EQ( lock2.mutex(), nullptr );
    BOOST_TEST( !lock2.owns_lock() );
    BOOST_TEST( !lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );

    BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
    BOOST_TEST_EQ( sp2.shared_lock_count(), 0 );
  }

  {
    // lhs bound, locked
    // rhs bound, locked

    dummy_lock sp, sp2;

    shared_lock_type lock( sp );
    shared_lock_type lock2( sp2 );

    lock2 = std::move( lock );
    BOOST_TEST_EQ( lock2.mutex(), std::addressof( sp ) );
    BOOST_TEST( lock2.owns_lock() );
    BOOST_TEST( lock2 );

    BOOST_TEST_EQ( lock.mutex(), nullptr );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
    BOOST_TEST_EQ( sp2.shared_lock_count(), 1 );
  }

  {
    // self-assign, locked

    dummy_lock sp;
    shared_lock_type lock( sp );

    lock = std::move( lock );
    BOOST_TEST_EQ( lock.mutex(), std::addressof( sp ) );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );

    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  }
}

void lock() {
  {
    dummy_lock sp;
    shared_lock_type lock( sp, std::defer_lock );

    lock.lock();
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );
    BOOST_TEST_THROWS( lock.lock(), std::system_error );

    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  }

  {
    dummy_lock sp;
    shared_lock_type lock;
    BOOST_TEST_THROWS( lock.lock(), std::system_error );
    BOOST_TEST_EQ( sp.shared_lock_count(), 0 );
  }
}

void unlock() {
  dummy_lock sp;

  shared_lock_type lock( sp );
  lock.unlock();

  BOOST_TEST( !lock.owns_lock() );
  BOOST_TEST( !lock );

  BOOST_TEST_EQ( sp.shared_lock_count(), 1 );

  sp.lock();
  sp.unlock();

  BOOST_TEST_THROWS( lock.unlock(), std::system_error );
}

void try_lock() {
  dummy_lock sp;

  {
    shared_lock_type lock( sp, std::defer_lock );
    BOOST_TEST( lock.try_lock() );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );

    BOOST_TEST( !sp.try_lock() );

    BOOST_TEST_THROWS( lock.try_lock(), std::system_error );

    lock.unlock();

    sp.lock();

    BOOST_TEST( !lock.try_lock() );
    BOOST_TEST( !lock.owns_lock() );
    BOOST_TEST( !lock );

    sp.unlock();
  }

  {
    shared_lock_type lock;
    BOOST_TEST_THROWS( lock.try_lock(), std::system_error );
  }
}

void swap() {
  dummy_lock sp, sp2;

  {
    shared_lock_type lock( sp );
    shared_lock_type lock2( sp2 );

    STATIC_ASSERT( noexcept( lock.swap( lock2 ) ) );
    STATIC_ASSERT( noexcept( swap( lock, lock2 ) ) );

    lock.swap( lock2 );

    BOOST_TEST_EQ( lock.mutex(), &sp2 );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );

    BOOST_TEST_EQ( lock2.mutex(), &sp );
    BOOST_TEST( lock2.owns_lock() );
    BOOST_TEST( lock2 );

    swap( lock, lock2 );

    BOOST_TEST_EQ( lock.mutex(), &sp );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );

    BOOST_TEST_EQ( lock2.mutex(), &sp2 );
    BOOST_TEST( lock2.owns_lock() );
    BOOST_TEST( lock2 );

    BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
    BOOST_TEST_EQ( sp2.shared_lock_count(), 1 );
  }

  {
    shared_lock_type lock( sp );
    lock.swap( lock );

    BOOST_TEST_EQ( lock.mutex(), &sp );
    BOOST_TEST( lock.owns_lock() );
    BOOST_TEST( lock );
  }
}

void release() {
  dummy_lock sp;

  shared_lock_type lock( sp );

  dummy_lock* pm = lock.release();

  BOOST_TEST_EQ( pm, &sp );
  BOOST_TEST_EQ( lock.mutex(), nullptr );
  BOOST_TEST( !lock.owns_lock() );
  BOOST_TEST( !lock );

  pm = lock.release();

  BOOST_TEST_EQ( lock.mutex(), nullptr );
  BOOST_TEST( !lock.owns_lock() );
  BOOST_TEST( !lock );

  BOOST_TEST_EQ( sp.shared_lock_count(), 1 );
  sp.unlock_shared();
}

} // namespace

int main() {
  sanity_tests();

  default_constructor();
  locking_construtor();
  deferred_constructor();
  try_lock_constructor();
  adopt_lock_constructor();

  move_constructor();
  move_assignment();

  lock();
  unlock();
  try_lock();

  swap();
  release();

  return boost::report_errors();
}

#ifdef BOOST_NO_EXCEPTIONS

namespace boost
{

BOOST_NORETURN void throw_exception( std::exception const& ex, boost::source_location const& loc )
{
    std::fprintf( stderr, "Exception '%s' at %s:%d\n", ex.what(), loc.file_name(), loc.line() );
    std::terminate();
}

} // namespace boost

#endif
