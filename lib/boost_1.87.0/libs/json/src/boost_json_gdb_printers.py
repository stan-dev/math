#
# Copyright (c) 2024 Dmitry Arkhipov (grisumbras@yandex.ru)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

import gdb
import gdb.printing


class PrettyPrinterBuilder():
    '''decorator that accumulates pretty printers for types'''

    def __init__(self, name=None):
        self.result = gdb.printing.RegexpCollectionPrettyPrinter(
            name or 'Boost.JSON')

    def __call__(self, pp=None, ns=None, template=False):
        if pp is None:
            def decorator(pp):
                return self(pp, ns=ns, template=template)
            return decorator

        typename = getattr(pp, '__name__')
        ns = ns or 'boost::json'
        self.result.add_printer(
            typename,
            '^{ns}::{typename}{marker}'.format(
                ns=ns,
                typename=typename,
                marker='<'if template else '$'),
            pp)
        return pp

pretty_printer = PrettyPrinterBuilder()


class static_property:
    '''decorator for lazy evaluation of static members'''

    def __init__(self, wrapped):
        self.wrapped = wrapped
        self.__name__ = wrapped.__name__
        self.__doc__ = wrapped.__doc__

    def __get__(self, inst, objtype=None):
        val = self.wrapped()
        if objtype is not None:
            setattr(objtype, self.wrapped.__name__, val)
        return val


@pretty_printer
class storage_ptr:
    '''boost::json::storage_ptr pretty printer'''

    @static_property
    def _size_t():
        return gdb.lookup_type('std::size_t')

    @static_property
    def _shared_resource_ptr():
        return gdb.lookup_type(
            'boost::json::detail::shared_resource').pointer()

    @static_property
    def _memory_resource_ptr():
        return gdb.lookup_type(
            'boost::container::pmr::memory_resource').pointer()

    def __init__(self, val):
        self.val = val['i_']

    def to_string(self):
        items = []

        if bool(self.val & 2):
            items.append('trivial')
        shared = bool(self.val & 1)
        if shared:
            items.append('shared')

        pointer = self.val & ~3
        if not pointer:
            items.append('resource=default')
        else:
            if shared:
                resource = pointer.cast(
                    self._shared_resource_ptr).dereference()
                impl_t = resource.dynamic_type
                impl_ptr = impl_t.pointer()
                resource = resource.address.cast(impl_ptr).dereference()

                items.append(
                    'refs=%s' % resource['refs'].cast(self._size_t))

                resource = resource['t']
            else:
                resource = pointer.cast(
                    self._memory_resource_ptr).dereference()


            derived_t = resource.dynamic_type
            derived_ptr = derived_t.pointer()
            resource = resource.address.cast(derived_ptr).dereference()
            items.append('resource=%s' % resource)

        return 'storage_ptr [' + ', '.join(items) + ']'


@pretty_printer
class monotonic_resource:
    '''boost::json::monotonic_resource pretty printer'''

    def __init__(self, val):
        self.val = val

    def to_string(self):
        buffer = self.val['buffer_']
        buffer = buffer['p'] - ( int(buffer['size']) - int(buffer['avail']) )

        head = self.val['head_'].dereference()

        block = head['p'] - ( int(head['size']) - int(head['avail']) )

        upstream = self.val['upstream_']

        items = []
        items.append('buffer=%s' % buffer)
        items.append('block=%s' % block)
        items.append('head=%s' % head['p'])
        items.append('free=%s' % head['avail'])
        if upstream['i_'] != 0:
            items.append('upstream=%s' % upstream)

        return 'monotonic_resource [%s]' % ', '.join(items)


@pretty_printer
class static_resource:
    '''boost::json::static_resource pretty printer'''

    def __init__(self, val):
        self.val = val

    def to_string(self):
        buf = self.val['p_'] - ( int(self.val['size_']) - int(self.val['n_']) )
        return 'static_resource [buffer=%s, head=%s, free=%s]' % (
            buf, self.val['p_'], self.val['n_'])


@pretty_printer
class string:
    '''boost::json::string pretty printer'''

    @static_property
    def _char_const_ptr():
        return gdb.lookup_type('char').const().pointer()

    def __init__(self, val):
        self.impl = val['impl_']

    def display_hint(self):
        return 'string'

    def to_string(self):
        kind = self.impl['s_']['k']
        if kind == self.impl['short_string_']:
            sbo_size = self.impl['sbo_chars_']
            size = sbo_size - self.impl['s_']['buf'][sbo_size]
            pointer = self.impl['s_']['buf']
        elif kind == self.impl['key_string_']:
            size = self.impl['k_']['n']
            pointer = self.impl['k_']['s']
        else:
            size = self.impl['p_']['t'].dereference()['size']
            pointer = self.impl['p_']['t']
            pointer += 1
            pointer = pointer.cast(self._char_const_ptr)
        return pointer.lazy_string(length=size)


@pretty_printer
class array:
    '''boost::json::array pretty printer'''

    @static_property
    def _value_ptr():
        return gdb.lookup_type('boost::json::value').pointer()

    def __init__(self, val):
        self.val = val

    def display_hint(self):
        return 'array'

    def to_string(self):
        capacity = int(self.val['t_'].dereference()['capacity'])
        return 'array [size={0}, capacity={1}]'.format(
            self.num_children(), capacity)

    def num_children(self):
        return int(self.val['t_'].dereference()['size'])

    def children(self):
        for i in range(0, self.num_children()):
            yield self.child(i)

    def child(self, n):
        table = (self.val['t_'] + 1).cast(self._value_ptr)
        return str(n), table[n]


@pretty_printer
class key_value_pair:
    '''boost::json::key_value_pair pretty printer'''

    @staticmethod
    def _pair(kv):
        return kv['key_'].lazy_string(length=kv['len_']), kv['value_']

    def __init__(self, val):
        self.val = val

    def to_string(self):
        k, v = self._pair(self.val)
        return '[%s] = %s' % (k.value(), v)


@pretty_printer
class object:
    '''boost::json::object pretty printer'''

    def __init__(self, val):
        self.val = val
        self.kv_ptr = gdb.lookup_type('boost::json::key_value_pair').pointer()

    def display_hint(self):
        return 'map'

    def to_string(self):
        capacity = int(self.val['t_'].dereference()['capacity'])
        return 'object [size={}, capacity={}]'.format(
            self.num_children(), capacity)

    def num_children(self):
        return int(self.val['t_'].dereference()['size'])

    def children(self):
        table = (self.val['t_'] + 1).cast(self.kv_ptr)
        for i in range(0, self.num_children()):
            k, v = key_value_pair._pair(table[i])
            yield str(2 * i), k
            yield str(2 * i + 1), v


@pretty_printer
class value:
    '''boost::json::value pretty printer'''

    @static_property
    def _kind_t():
        result = gdb.lookup_type('boost::json::kind')
        return dict(result.items())

    def __init__(self, val):
        self.val = val

    def to_string(self):
        kind = self.val['sca_']['k']

        if self._compare_kind(kind, 'null'):
            return 'null'

        elif self._compare_kind(kind, 'bool_'):
            return self.val['sca_']['b']

        elif self._compare_kind(kind, 'int64'):
            return self.val['sca_']['i']

        elif self._compare_kind(kind, 'uint64'):
            return self.val['sca_']['u']

        elif self._compare_kind(kind, 'double_'):
            return self.val['sca_']['d']

        elif self._compare_kind(kind, 'array'):
            return self.val['arr_']

        elif self._compare_kind(kind, 'object'):
            return self.val['obj_']

        else:
            return self.val['str_']

    def _compare_kind(self, kind, name):
        return kind == self._kind_t['boost::json::kind::' + name].enumval


def register(obj_file):
    mod = obj_file or gdb
    for printer in getattr(mod, 'pretty_printers', []):
        if getattr(printer, 'name') == pretty_printer.result.name:
            return

    gdb.printing.register_pretty_printer(
        obj_file,
        pretty_printer.result)


if __name__ == '__main__':
    register(gdb.current_objfile())
