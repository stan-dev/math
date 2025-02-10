/*
 * Copyright 1993, 1995 Christopher Seiwald.
 *
 * This file is part of Jam - see jam.c for Copyright information.
 */

/*  This file is ALSO:
 *  Copyright 2022 René Ferdinand Rivera Morell
 *  Copyright 2001-2004 David Abrahams.
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE.txt or
 * https://www.bfgroup.xyz/b2/LICENSE.txt)
 */

/*
 * lists.h - the LIST structure and routines to manipulate them
 *
 * The whole of jam relies on lists of objects as a datatype. This module, in
 * conjunction with object.c, handles these relatively efficiently.
 *
 * Structures defined:
 *
 *  LIST - list of OBJECTs
 *  LOL - list of LISTs
 *
 * External routines:
 *
 *  list_append() - append a list onto another one, returning total
 *  list_new() - tack an object onto the end of a list of objects
 *  list_copy() - copy a whole list of objects
 *  list_sublist() - copy a subset of a list of objects
 *  list_free() - free a list of objects
 *  list_print() - print a list of objects to stdout
 *  list_length() - return the number of items in the list
 *
 *  lol_init() - initialize a LOL (list of lists)
 *  lol_add() - append a LIST onto an LOL
 *  lol_free() - free the LOL and its LISTs
 *  lol_get() - return one of the LISTs in the LOL
 *  lol_print() - debug print LISTS separated by ":"
 */

#ifndef B2_LISTS_H
#define B2_LISTS_H

#include "config.h"

#include "jam.h"
#include "object.h"
#include "output.h"
#include "value.h"

/*
 * LIST - list of strings
 */

struct LIST
{
	union
	{
		int32_t size;
		struct LIST * next;
		OBJECT * align;
	} impl;

	LIST() { this->impl.next = nullptr; }
};

typedef LIST * list_ptr;
typedef OBJECT ** LISTITER;

/*
 * LOL - list of LISTs
 */

#define LOL_MAX 19
typedef struct _lol
{
	int32_t count;
	LIST * list[LOL_MAX];
} LOL;

LIST * list_new(OBJECT * value);
LIST * list_append(LIST * destination, LIST * source);
LIST * list_copy(LIST *);
LIST * list_copy_range(LIST * destination, LISTITER first, LISTITER last);
void list_free(LIST * head);
LIST * list_push_back(LIST * head, OBJECT * value);
void list_print(LIST *);
inline int32_t list_length(LIST * l) { return l ? l->impl.size : 0; }
LIST * list_sublist(LIST *, int32_t start, int32_t count);
LIST * list_pop_front(LIST *);
LIST * list_sort(LIST *);
LIST * list_unique(LIST * sorted_list);
int32_t list_in(LIST *, OBJECT * value);
LIST * list_reverse(LIST *);
int32_t list_cmp(LIST * lhs, LIST * rhs);
int32_t list_is_sublist(LIST * sub, LIST * l);
void list_done();

#define L0 ((LIST *)0)

inline LISTITER list_begin(LIST * l)
{
	return l ? (LISTITER)((char *)l + sizeof(LIST)) : 0;
}
inline LISTITER list_end(LIST * l)
{
	return l ? list_begin(l) + l->impl.size : 0;
}
inline LISTITER list_next(LISTITER it) { return ((it) + 1); }
inline LISTITER list_prev(LISTITER it) { return ((it)-1); }
inline OBJECT *& list_item(LISTITER it) { return (*(it)); }
inline bool list_empty(LIST * l) { return ((l) == L0); }
inline OBJECT *& list_front(LIST * l) { return list_item(list_begin(l)); }
inline OBJECT *& list_item(LIST * l, int i)
{
	return list_begin(l)[i < 0 ? l->impl.size + i : i];
}
void lol_add_err();
inline void lol_add(LOL * lol, LIST * l)
{
	if (lol->count < LOL_MAX)
	{
		lol->list[lol->count++] = l;
		return;
	}
	lol_add_err();
}
inline void lol_init(LOL * lol) { lol->count = 0; }
inline void lol_free(LOL * lol)
{
	int32_t i;
	for (i = 0; i < lol->count; ++i) list_free(lol->list[i]);
	lol->count = 0;
}
inline LIST * lol_get(LOL * lol, int32_t i)
{
	return i < lol->count ? lol->list[i] : L0;
}
inline void lol_print(LOL * lol)
{
	int32_t i;
	for (i = 0; i < lol->count; ++i)
	{
		if (i) out_printf(" : ");
		list_print(lol->list[i]);
	}
}
void lol_build(LOL * lol, char const ** elements);

namespace b2 {

struct list_ref;

/* tag::reference[]

= `b2::list_cref`

Container of b2 values, that is a non-owning reference to the `LIST`. Mostly
follows random access container behavior.

== `b2::list_cref` Overview

[source,cpp]
----
end::reference[] */
// tag::reference[]
struct list_cref
{
	// types
	struct iterator;
	using size_type = int32_t;
	using value_type = OBJECT *;

	// construct/copy/destroy
	list_cref() = default;
	list_cref(const list_cref &) = default;
	list_cref(list_cref && other);
	explicit list_cref(LIST * l);
	list_cref & operator=(const list_cref &) = default;

	// iterators
	iterator begin() const;
	iterator end() const;

	// capacity
	bool empty() const B2_NOEXCEPT;
	size_type length() const B2_NOEXCEPT;
	size_type size() const B2_NOEXCEPT;

	// element access
	value_type & operator[](size_type i) const;

	// list operations
	bool contains(value_ref a) const;
	list_ref slice(size_type i, size_type j = -1) const;
	bool operator==(const list_cref & b) const;
	bool operator==(const list_ref & b) const;

	// data access
	LIST * data() const B2_NOEXCEPT;
	LIST * operator*() const B2_NOEXCEPT;

	protected:
	friend struct iterator;
	LIST * list_obj = nullptr;
};
// end::reference[]
/* tag::reference[]
----
end::reference[] */

struct list_cref::iterator
{
	using iterator_category = std::random_access_iterator_tag;
	using value_type = OBJECT *;
	using difference_type = std::ptrdiff_t;
	using pointer = value_type *;
	using reference = value_type &;

	inline explicit iterator(LISTITER i)
		: list_i(i)
	{}

	inline iterator operator++()
	{
		list_i = list_next(list_i);
		return *this;
	}
	inline iterator operator++(int)
	{
		iterator result { *this };
		list_i = list_next(list_i);
		return result;
	}
	inline iterator operator--()
	{
		list_i = list_prev(list_i);
		return *this;
	}
	inline iterator operator--(int)
	{
		iterator result { *this };
		list_i = list_prev(list_i);
		return result;
	}
	inline bool operator==(iterator other) const
	{
		return list_i == other.list_i;
	}
	inline bool operator!=(iterator other) const
	{
		return list_i != other.list_i;
	}
	inline reference operator*() const { return list_item(list_i); }
	inline pointer operator->() const { return &list_item(list_i); }
	inline difference_type operator-(const iterator & b) const
	{
		return list_i - b.list_i;
	}
	inline iterator operator+(difference_type count) const
	{
		return iterator(list_i + count);
	}
	inline iterator operator-(difference_type count) const
	{
		return iterator(list_i - count);
	}

	private:
	LISTITER list_i;

	// MSVC secure/debug std asks for comparisons of iterators.
	friend inline bool operator<(const iterator & a, const iterator & b)
	{
		return a.list_i < b.list_i;
	}
};

/* tag::reference[]
== `b2::list_cref` Construct/Copy/Destroy
end::reference[] */

/* tag::reference[]
=== `b2::list_cref::list_cref`

[source,cpp]
----
inline list_cref::list_cref(list_cref && other)
inline list_cref::list_cref(LIST * l)
----
end::reference[] */
inline list_cref::list_cref(list_cref && other)
	: list_obj(std::move(other.list_obj))
{}
inline list_cref::list_cref(LIST * l)
	: list_obj(l)
{}

/* tag::reference[]
== `b2::list_cref` Iterators
end::reference[] */

/* tag::reference[]
=== `b2::list_cref::begin`

[source,cpp]
----
inline list_cref::iterator list_cref::begin() const
----
end::reference[] */
inline list_cref::iterator list_cref::begin() const
{
	return iterator(list_begin(list_obj));
}

/* tag::reference[]
=== `b2::list_cref::end`

[source,cpp]
----
inline list_cref::iterator list_cref::end() const
----
end::reference[] */
inline list_cref::iterator list_cref::end() const
{
	return iterator(list_end(list_obj));
}

/* tag::reference[]
== `b2::list_cref` Capacity
end::reference[] */

/* tag::reference[]
=== `b2::list_cref::empty`

[source,cpp]
----
inline bool list_cref::empty() const B2_NOEXCEPT
----
end::reference[] */
inline bool list_cref::empty() const B2_NOEXCEPT
{
	return list_empty(list_obj) || length() == 0;
}

/* tag::reference[]
=== `b2::list_cref::length`

[source,cpp]
----
inline list_cref::size_type list_cref::length() const B2_NOEXCEPT
inline list_cref::size_type list_cref::size() const B2_NOEXCEPT
----
end::reference[] */
inline list_cref::size_type list_cref::length() const B2_NOEXCEPT
{
	return list_length(list_obj);
}
inline list_cref::size_type list_cref::size() const B2_NOEXCEPT
{
	return length();
}

/* tag::reference[]
== `b2::list_cref` Element Access
end::reference[] */

/* tag::reference[]
=== `b2::list_cref::operator[]`

[source,cpp]
----
inline list_cref::value_type & list_cref::operator[](
	list_cref::size_type i) const
----
end::reference[] */
inline list_cref::value_type & list_cref::operator[](
	list_cref::size_type i) const
{
	return list_item(list_obj, i);
}

/* tag::reference[]
== `b2::list_cref` List Operations
end::reference[] */

/* tag::reference[]
=== `b2::list_cref::contains`

[source,cpp]
----
inline bool list_cref::contains(value_ref a) const
----
end::reference[] */
inline bool list_cref::contains(value_ref a) const
{
	for (auto b : *this)
		if (a == b) return true;
	return false;
}

/* tag::reference[]
=== `b2::list_cref::slice`

[source,cpp]
----
inline list_ref list_cref::slice(
	list_cref::size_type i, list_cref::size_type j) const
----
end::reference[] */

/* tag::reference[]
=== `b2::list_cref::operator==`

[source,cpp]
----
inline bool list_cref::operator==(const list_cref & b) const
inline bool list_cref::operator==(const list_ref & b) const
----
end::reference[] */
inline bool list_cref::operator==(const list_cref & b) const
{
	if (length() != b.length()) return false;
	iterator b_i = b.begin();
	for (value_ref v : *this)
		if (v != value_ref(*(b_i++))) return false;
	return true;
}

/* tag::reference[]
== `b2::list_cref` Data Access
end::reference[] */

/* tag::reference[]
=== `b2::list_cref::operator==`

[source,cpp]
----
inline LIST * list_cref::data() const B2_NOEXCEPT
inline LIST * list_cref::operator*() const B2_NOEXCEPT
----
end::reference[] */
inline LIST * list_cref::data() const B2_NOEXCEPT { return list_obj; }
inline LIST * list_cref::operator*() const B2_NOEXCEPT { return data(); }

/* tag::reference[]

= `b2::list_ref`

Container of b2 values, that is an owning reference to the `LIST`. Mostly
follows random access container behavior. And as an owning reference will
allocate, copy, move `LIST` objects as needed.

== `b2::list_ref` Overview

[source,cpp]
----
end::reference[] */
// tag::reference[]
struct list_ref : private list_cref
{
	// types
	using list_cref::iterator;
	using list_cref::size_type;
	using list_cref::value_type;

	using list_cref::begin;
	using list_cref::end;
	using list_cref::empty;
	using list_cref::length;
	using list_cref::size;
	using list_cref::operator[];
	using list_cref::contains;
	using list_cref::operator==;
	using list_cref::data;
	using list_cref::operator*;

	// construct/copy/destroy
	list_ref() = default;
	list_ref(list_ref && other);
	list_ref(const list_cref & other);
	list_ref(const list_ref & other);
	explicit list_ref(value_ref o);
	explicit list_ref(LIST * l, bool own = false);
	list_ref(iterator i, const iterator & e);
	~list_ref();

	// modifiers
	LIST * release();
	void reset(LIST * new_list = nullptr);
	list_ref & append(const list_ref & other);
	list_ref & append(list_cref other);
	list_ref & operator+(const list_ref & other);
	list_ref & operator+(const list_cref & other);
	list_ref & push_back(OBJECT * value);
	template <typename... T>
	list_ref & push_back(T... value);
	template <typename T>
	list_ref & operator+(T value);
	list_ref & pop_front();
	list_ref & operator=(list_ref && other);

	// list operations
	inline list_ref & slice(size_type i, size_type j = -1);
	inline list_cref cref() const;
};
// end::reference[]
/* tag::reference[]
----
end::reference[] */

inline list_ref list_cref::slice(
	list_cref::size_type i, list_cref::size_type j) const
{
	return list_ref(
		list_sublist(list_obj, i, (j < 0 ? length() + j : j) - i + 1), true);
}

inline bool list_cref::operator==(const list_ref & b) const
{
	return *this == list_cref(*b);
}

/* tag::reference[]
== `b2::list_ref` Construct/Copy/Destroy
end::reference[] */

/* tag::reference[]
=== `b2::list_ref::list_cref`

[source,cpp]
----
inline list_ref::list_ref(list_ref && other) // <1>
inline list_ref::list_ref(const list_cref & other) // <2>
inline list_ref::list_ref(const list_ref & other) // <2>
inline list_ref::list_ref(value_ref o) // <3>
inline list_ref::list_ref(LIST * l, bool own) // <4>
inline list_ref::list_ref(iterator i, const iterator & e) // <5>
----
<1> The data for the list is moved from `other`.
<2> Makes a copy of the list from `other`.
<3> Makes a new list with the initial given element value.
<4> If `own == true` takes ownership of the data `l`, otherwise makes a copy of
	it.
<5> Fills the new list with the elements from the `[i,e)` range.
end::reference[] */
inline list_ref::list_ref(list_ref && other)
	: list_cref(other.release())
{}
inline list_ref::list_ref(const list_cref & other)
	: list_cref(list_copy(*other))
{}
inline list_ref::list_ref(const list_ref & other)
	: list_cref(list_copy(other.list_obj))
{}
inline list_ref::list_ref(value_ref o)
	: list_cref(list_new(o))
{}
inline list_ref::list_ref(LIST * l, bool own)
	: list_cref(own ? l : list_copy(l))
{}
inline list_ref::list_ref(iterator i, const iterator & e)
{
	for (; i != e; ++i) this->push_back(value::copy(*i));
}
inline list_ref::~list_ref() { reset(); }

/* tag::reference[]
== `b2::list_ref` Modifiers
end::reference[] */

/* tag::reference[]
=== `b2::list_ref::release`

[source,cpp]
----
inline LIST * list_ref::release()
----

Returns the list data relinquishing ownership of it. This list is left in an
empty valid state.
end::reference[] */
inline LIST * list_ref::release()
{
	LIST * r = list_obj;
	list_obj = nullptr;
	return r;
}

/* tag::reference[]
=== `b2::list_ref::reset`

[source,cpp]
----
inline void list_ref::reset(LIST * new_list)
----

Replaces the list data with the given `new_list`. The current list is freed
along with the elements in the list.
end::reference[] */
inline void list_ref::reset(LIST * new_list)
{
	std::swap(list_obj, new_list);
	if (new_list) list_free(new_list);
}

/* tag::reference[]
=== `b2::list_ref::append`

[source,cpp]
----
inline list_ref & list_ref::append(const list_ref & other)
inline list_ref & list_ref::append(list_cref other)
inline list_ref & list_ref::operator+(const list_ref & other)
inline list_ref & list_ref::operator+(const list_cref & other)
----

Adds the elements from the `other` list, making copies, to the end of this list.
All the functions return a reference to the list to allow for chaining. For
example: `list_ref() + "one" + "two"`.
end::reference[] */
inline list_ref & list_ref::append(const list_ref & other)
{
	list_obj = list_append(list_obj, list_copy(*other));
	return *this;
}
inline list_ref & list_ref::append(list_cref other)
{
	list_obj = list_append(list_obj, list_copy(*other));
	return *this;
}
inline list_ref & list_ref::operator+(const list_ref & other)
{
	return append(other);
}
inline list_ref & list_ref::operator+(const list_cref & other)
{
	return append(other);
}

/* tag::reference[]
=== `b2::list_ref::push_back`

[source,cpp]
----
inline list_ref & list_ref::push_back(OBJECT * value)
template <typename... T>
inline list_ref & list_ref::push_back(T... value) // <1>
template <typename... T>
inline list_ref & list_ref::operator+(T value) // <2>
----
<1> Adds a value constructed from the given arguments. I.e. by calling
	`value::make(value...)`.
<2> Adds the `value` that is convertible to `value_type`.

Adds a single value to the end of the list. The list is returned to allow for
chaining.
end::reference[] */
inline list_ref & list_ref::push_back(OBJECT * value)
{
	list_obj = list_push_back(list_obj, value);
	return *this;
}
template <typename... T>
inline list_ref & list_ref::push_back(T... value)
{
	list_obj = list_push_back(list_obj, value::make(value...));
	return *this;
}
template <typename T>
inline list_ref & list_ref::operator+(T value)
{
	return push_back(value);
}

/* tag::reference[]
=== `b2::list_ref::pop_front`

[source,cpp]
----
inline list_ref & list_ref::pop_front()
----

end::reference[] */
inline list_ref & list_ref::pop_front()
{
	list_obj = list_pop_front(list_obj);
	return *this;
}

/* tag::reference[]
=== `b2::list_ref::operator=`, assign

[source,cpp]
----
inline list_ref & list_ref::operator=(list_ref && other) // <1>
----
<1> Moves the data from `other` to this list.

end::reference[] */
inline list_ref & list_ref::operator=(list_ref && other)
{
	reset(other.list_obj);
	other.list_obj = nullptr;
	return *this;
}

/* tag::reference[]
=== `b2::list_ref::slice`

[source,cpp]
----
inline list_ref & list_ref::slice(size_type i, size_type j)
----

Replaces, in-place, the list with the indicated subrange slice `[i,j]`. Negative
values of `j` index from the `end()` position.

end::reference[] */
inline list_ref & list_ref::slice(size_type i, size_type j)
{
	list_obj = list_sublist(list_obj, i, (j < 0 ? length() + j : j) - i + 1);
	return *this;
}

inline list_cref list_ref::cref() const { return list_cref(list_obj); }

/* tag::reference[]

= `b2::lists`

Container of a "list of list" that is owning an instance of `LOL` object. The
interface allows for inline composition of the `LOL`.

== `b2::lists` Overview

[source,cpp]
----
end::reference[] */
// tag::reference[]
struct lists
{
	// types
	using size_type = int32_t;

	// construct/copy/destroy
	lists();
	lists(lists && other);
	~lists();

	// capacity
	bool empty() const B2_NOEXCEPT;
	size_type size() const B2_NOEXCEPT;
	size_type length() const B2_NOEXCEPT;
	size_type max_size() const B2_NOEXCEPT;
	size_type capacity() const B2_NOEXCEPT;

	// element access
	list_cref operator[](size_type i) const;

	// modifiers
	void push_back(const list_cref & l);
	void push_back(list_ref && l);
	lists & operator|(const list_cref & l);
	lists & operator|(LIST * l);
	lists & operator|(list_ref && l);
	lists & operator|=(list_ref && l);
	lists & append(const lists & lol);
	lists & operator|(const lists & lol);
	void swap(LOL & other);
	void clear();

	// display
	void print() const;

	// data access
	LOL * data() const B2_NOEXCEPT;
	operator LOL *() const B2_NOEXCEPT;

	private:
	mutable LOL lol;
};
// end::reference[]
/* tag::reference[]
----
end::reference[] */

/* tag::reference[]
== `b2::lists` Construct/Copy/Destroy
end::reference[] */

/* tag::reference[]
=== `b2::lists::lists`

[source,cpp]
----
inline lists::lists()
inline lists::lists(lists && other)
----
end::reference[] */
inline lists::lists() { lol_init(&lol); }
inline lists::lists(lists && other)
{
	lol_init(&lol);
	std::swap(lol, other.lol);
}

/* tag::reference[]
=== `b2::lists::~lists`

[source,cpp]
----
inline lists::~lists()
----
end::reference[] */
inline lists::~lists() { lol_free(&lol); }

/* tag::reference[]
== `b2::lists` Capacity
end::reference[] */
inline bool lists::empty() const B2_NOEXCEPT { return length() == 0; }
inline lists::size_type lists::size() const B2_NOEXCEPT { return lol.count; }
inline lists::size_type lists::length() const B2_NOEXCEPT { return lol.count; }
inline lists::size_type lists::max_size() const B2_NOEXCEPT { return LOL_MAX; }
inline lists::size_type lists::capacity() const B2_NOEXCEPT { return LOL_MAX; }

/* tag::reference[]
== `b2::lists` Element Access
end::reference[] */

/* tag::reference[]
=== `b2::lists::operator[]`

[source,cpp]
----
inline list_cref lists::operator[](int32_t i) const
----

Returns a constant reference, i.e. `list_cref`, of the list at the given `i`
index.
end::reference[] */
inline list_cref lists::operator[](list_cref::size_type i) const
{
	return list_cref(lol_get(&lol, i));
}

/* tag::reference[]
== `b2::lists` Modifiers
end::reference[] */

/* tag::reference[]
=== `b2::lists::push_back`

[source,cpp]
----
inline void lists::push_back(const list_cref & l)
inline void lists::push_back(list_ref && l)
----

Adds the given list to the end of the `LOL`.
end::reference[] */
inline void lists::push_back(const list_cref & l)
{
	lol_add(&lol, list_copy(*l));
}
inline void lists::push_back(list_ref && l)
{
	list_ref to_add(std::move(l));
	lol_add(&lol, to_add.release());
}

/* tag::reference[]
=== `b2::lists::operator|`

[source,cpp]
----
inline lists & lists::operator|(const list_cref & l)
inline lists & lists::operator|(LIST * l)
inline lists & lists::operator|(list_ref && l)
inline lists & lists::operator|=(list_ref && l)
----

Adds the given list to the end of the `LOL`. This returns this object making it
possible to chain additions into a single statement.
end::reference[] */
inline lists & lists::operator|(const list_cref & l)
{
	push_back(l);
	return *this;
}
inline lists & lists::operator|(LIST * l)
{
	push_back(list_cref(l));
	return *this;
}
inline lists & lists::operator|(list_ref && l)
{
	push_back(std::move(l));
	return *this;
}
inline lists & lists::operator|=(list_ref && l)
{
	push_back(std::move(l));
	return *this;
}

inline lists & lists::append(const lists & lol)
{
	for (size_type i = 0; i < lol.length(); ++i) push_back(lol[i]);
	return *this;
}
inline lists & lists::operator|(const lists & lol) { return append(lol); }

/* tag::reference[]
=== `b2::lists::swap`

[source,cpp]
----
inline void lists::swap(LOL & other)
----

Swaps the data from `other` with this.
end::reference[] */
inline void lists::swap(LOL & other) { std::swap(lol, other); }

/* tag::reference[]
=== `b2::lists::clear`

[source,cpp]
----
inline void lists::clear()
----

Deallocate any list items and resets the length to zero.
end::reference[] */
inline void lists::clear() { lol_free(&lol); }

/* tag::reference[]
== `b2::lists` Element Access
end::reference[] */

/* tag::reference[]
=== `b2::lists::print`

[source,cpp]
----
inline void lists::print() const
----

Outputs, to cout, the lists as colon (`:`) separated and space separated
elements.
end::reference[] */
inline void lists::print() const { lol_print(&lol); }

/* tag::reference[]
== `b2::lists` Data Access
end::reference[] */

/* tag::reference[]
=== `b2::lists::data`

[source,cpp]
----
inline LOL * lists::data() const B2_NOEXCEPT
inline lists::operator LOL *() const B2_NOEXCEPT
----

Returns the underlying `LOL` object.
end::reference[] */
inline LOL * lists::data() const B2_NOEXCEPT { return &lol; }
inline lists::operator LOL *() const B2_NOEXCEPT { return data(); }

} // namespace b2

#endif
