#
# Copyright (c) 2023 alandefreitas (alandefreitas@gmail.com)
#
# Distributed under the Boost Software License, Version 1.0.
# https://www.boost.org/LICENSE_1_0.txt
#

import gdb
import urllib.parse


class utils:
    @staticmethod
    def resolve_type(t):
        if t.code == gdb.TYPE_CODE_REF:
            t = t.target()
        t = t.unqualified().strip_typedefs()
        typename = t.tag
        if typename is None:
            return None
        return t

    @staticmethod
    def resolved_typename(val):
        t = val.type
        t = utils.resolve_type(t)
        if t is not None:
            return str(t)
        else:
            return str(val.type)

    @staticmethod
    def pct_decode(s: str):
        return urllib.parse.unquote(s)

    sv_pool = []

    @staticmethod
    def make_string_view(cstr: gdb.Value, n: int):
        sv_ptr: gdb.Value = gdb.parse_and_eval('malloc(sizeof(class boost::core::basic_string_view<char>))')
        sv_ptr_str = cstr.format_string(format='x')
        gdb.execute(
            f'call ((boost::core::basic_string_view<char>*){sv_ptr})->basic_string_view((const char*){sv_ptr_str}, {n})',
            to_string=True)
        sv = gdb.parse_and_eval(f'*((boost::core::basic_string_view<char>*){sv_ptr})')
        copy: gdb.Value = sv
        utils.sv_pool.append(sv_ptr)
        if len(utils.sv_pool) > 5000:
            gdb.execute(f'call free({utils.sv_pool[0]})', to_string=True)
            utils.sv_pool.pop(0)
        return copy

    pct_sv_pool = []

    @staticmethod
    def make_pct_string_view(cstr: gdb.Value, n: int):
        sv_ptr: gdb.Value = gdb.parse_and_eval('malloc(sizeof(class boost::urls::pct_string_view))')
        sv_ptr_str = cstr.format_string(format='x')
        gdb.execute(
            f'call ((boost::urls::pct_string_view*){sv_ptr})->pct_string_view((const char*){sv_ptr_str}, {n})',
            to_string=True)
        sv = gdb.parse_and_eval(f'*((boost::urls::pct_string_view*){sv_ptr})')
        copy: gdb.Value = sv
        utils.pct_sv_pool.append(sv_ptr)
        if len(utils.sv_pool) > 5000:
            gdb.execute(f'call free({utils.pct_sv_pool[0]})', to_string=True)
            utils.sv_pool.pop(0)
        return copy


class StringViewPrinter:
    def __init__(self, value):
        self.value = value

    def to_string(self):
        s = str(self.value['p_'])
        ptr, s = s.split(' ', 1)
        s = s[:self.value['n_'] + 1] + '"'
        return f"{ptr} {s} ({self.value['n_']})"
        # return f"{{ptr} {s} ({self.value['n_']})"


class PctStringViewPrinter:
    def __init__(self, value):
        self.value = value

    def to_string(self):
        if self.value['dn_'] != self.value['s_']['n_']:
            sd = utils.pct_decode(str(self.value['s_']['p_']))
            ptr, sd = sd.split(' ', 1)
            sd = sd[:self.value['dn_'] + 1] + '"'
            return f"{self.value['s_']}: {sd} ({self.value['dn_']})"
        else:
            return f"{self.value['s_']}"


class UrlImplPrinter:
    parts = ['scheme', 'user', 'pass', 'host', 'port', 'path', 'query', 'fragment']
    grammar_part = ['scheme ":"', '"//" user', '":" pass "@"', 'host', '":" port', 'path', '"?" query',
                    '"#" fragment']
    part_delimiters = [1, 2, 2, 0, 1, 0, 1, 1]
    show_offset_parts = False

    def __init__(self, value, prefix=None, print_cstr=None):
        self.value = value
        self.prefix = '' if prefix is None else prefix + '::'
        self.print_cstr = True if print_cstr is None else print_cstr

    def to_string(self):
        return self.value['cs_']

    def children(self):
        if self.print_cstr:
            yield 'buffer', self.value['cs_']

        yield 'size', self.value['offset_'][7]

        # Parts
        part_name = self.grammar_part[0]
        component_name = self.parts[0]
        part_size = self.value['offset_'][0]
        if part_size != 0:
            if self.show_offset_parts:
                yield f"{part_name}", utils.make_string_view(self.value['cs_'], part_size)
            else:
                yield f"{component_name}", utils.make_string_view(self.value['cs_'], part_size - 1)
                yield f"scheme id", self.value['scheme_']
        for i in range(7):
            cs_part = self.value['cs_'] + self.value['offset_'][i]
            part_name = self.grammar_part[i + 1]
            component_name = self.parts[i + 1]
            part_offset = self.value['offset_'][i]
            next_part_offset = self.value['offset_'][i + 1]
            part_size = next_part_offset - part_offset
            if self.show_offset_parts:
                yield f"{part_name}", utils.make_pct_string_view(cs_part, part_size)
            elif part_size != 0:
                component_part_offset = {
                    'user': 2,
                    'pass': 1,
                    'port': 1,
                    'query': 1,
                    'fragment': 1
                }
                component_part_size_diff = {
                    'user': 2,
                    'pass': 2,
                    'port': 1,
                    'query': 1,
                    'fragment': 1
                }
                off = component_part_offset[component_name] if component_name in component_part_offset else 0
                diff = component_part_size_diff[component_name] if component_name in component_part_size_diff else 0
                if part_size - diff > 0:
                    yield f"{component_name}", utils.make_pct_string_view(cs_part + off, part_size - diff)
                    if component_name == 'host':
                        yield f"host type", self.value['host_type_']
                        if not str(self.value['host_type_']).endswith('name'):
                            yield f"ip address", self.value['ip_addr_']
                    elif component_name == 'port':
                        yield f"port number", self.value['port_number_']
                    elif component_name == 'path':
                        yield f"# segments", self.value['nseg_']
                    elif component_name == 'query':
                        yield f"# params", self.value['nparam_']

        already_yielded = ['zero_', 'cs_', 'offset_', 'decoded_', 'nseg_', 'nparam_', 'scheme_',
                           'port_number_', 'ip_addr_', 'host_type_']
        for field in self.value.type.fields():
            if field.is_base_class:
                if not field.name.endswith('parts_base'):
                    yield self.prefix + field.name, self.value.cast(field.type)
            elif field.artificial:
                continue
            elif field.name not in already_yielded:
                name = field.name if not field.name.endswith('_') else field.name[:-1]
                yield name, self.value[field.name]


class UrlViewBasePrinter:
    def __init__(self, value, prefix=None, print_cstr=None):
        self.value = value
        self.prefix = '' if prefix is None else prefix + '::'
        self.print_cstr = True if print_cstr is None else print_cstr

    def to_string(self):
        return self.value['s_']

    def children(self):
        obj_type = self.value.type

        inline_impl: bool = self.value['pi_'] == self.value['impl_'].address
        if inline_impl:
            yield from UrlImplPrinter(self.value['impl_'], 'url_impl', print_cstr=self.print_cstr).children()
        else:
            yield 'reference', self.value['pi_']

        for field in obj_type.fields():
            if field.is_base_class:
                if not field.name.endswith('parts_base'):
                    yield self.prefix + field.name, self.value.cast(field.type)
            elif field.artificial:
                continue
            elif field.name not in ['pi_', 'impl_']:
                yield self.prefix + field.name, self.value[field.name]


class UrlBasePrinter:
    def __init__(self, value, prefix=None):
        self.value = value
        self.prefix = '' if prefix is None else prefix + '::'

    def to_string(self):
        return self.value['s_']

    def children(self):
        yield 'buffer', self.value['s_']
        yield 'capacity', self.value['cap_']
        for field in self.value.type.fields():
            if field.is_base_class or field.artificial:
                continue
            elif field.name not in ['s_', 'cap_']:
                yield self.prefix + field.name, self.value[field.name]
        for field in self.value.type.fields():
            if field.is_base_class:
                if field.name == 'boost::urls::url_view_base':
                    yield from UrlViewBasePrinter(self.value.cast(field.type), 'url_view_base',
                                                  print_cstr=False).children()
                else:
                    yield self.prefix + field.name, self.value.cast(field.type)


class UrlPrinter:
    def __init__(self, value):
        self.value = value

    def to_string(self):
        return self.value['pi_']['s_']

    def children(self):
        obj_type = self.value.type

        for field in obj_type.fields():
            if not field.is_base_class and not field.artificial:
                yield field.name, self.value[field.name]

        for field in obj_type.fields():
            if field.is_base_class:
                if field.name == 'boost::urls::url_base':
                    yield from UrlBasePrinter(self.value.cast(field.type), 'url_base').children()
                else:
                    yield field.name, self.value.cast(field.type)


class UrlViewPrinter:
    def __init__(self, value):
        self.value = value

    def to_string(self):
        return self.value['pi_']['cs_']

    def children(self):
        obj_type = self.value.type

        for field in obj_type.fields():
            if not field.is_base_class and not field.artificial:
                yield field.name, self.value[field.name]

        for field in obj_type.fields():
            if field.is_base_class:
                if field.name == 'boost::urls::url_view_base':
                    yield from UrlViewBasePrinter(self.value.cast(field.type), 'url_view_base',
                                                  print_cstr=True).children()
                else:
                    yield field.name, self.value.cast(field.type)


class Variant2Printer:
    def __init__(self, value):
        self.value = value

    def children(self):
        cur_value = self.value
        while cur_value is not None and not self.has_member(cur_value, 'st_') and not self.has_member(cur_value, 'ix_'):
            cur_value = self.get_impl(cur_value)

        if cur_value is not None:
            # Print the value of the st_ field
            index = cur_value['ix_']
            storage = cur_value['st_']
            i = 0
            storage_value = storage['first_']
            while i < index:
                storage = storage['rest_']
                storage_value = storage['first_']
                i += 1
            yield 'storage', storage_value
            yield 'index', index
        else:
            # We could not find the st_ field, so we just print the value
            obj_type = self.value.type
            for field in obj_type.fields():
                if not field.is_base_class and not field.artificial:
                    yield field.name, self.value[field.name]

            for field in obj_type.fields():
                if field.is_base_class:
                    yield field.name, self.value.cast(field.type)

    @staticmethod
    def has_member(value, member_name):
        obj_type = value.type
        for field in obj_type.fields():
            if field.name == member_name:
                return True
        return False

    # static function to get the value casted to the first base class type whose the substring before the first < ends with _impl
    @staticmethod
    def get_impl(value):
        obj_type = value.type
        for field in obj_type.fields():
            if field.is_base_class:
                if field.name.startswith('boost::variant2::detail::variant_'):
                    # Get substring of field.name before the first <
                    name = field.name.split('<', 1)[0]
                    if name.endswith('_impl'):
                        return value.cast(field.type)
        return None


class SystemResultPrinter:
    def __init__(self, value):
        self.value = value

    def children(self):
        yield from Variant2Printer(self.value['v_']).children()


class EnumPrinter:
    def __init__(self, value):
        self.value = value

    def to_string(self):
        s: str = self.value.format_string(raw=True)
        return s.rsplit(':', 1)[-1]


class UrlHostTypePrinter:
    def __init__(self, value):
        self.value = value

    def to_string(self):
        s: str = self.value.format_string(raw=True)
        s = s.rsplit(':', 1)[-1]
        return 'reg-name' if s == 'name' else s


if __name__ != "__main__":
    def lookup_function(val: gdb.Value):
        typename: str = utils.resolved_typename(val)

        if typename == 'boost::urls::url_view':
            return UrlViewPrinter(val)
        elif typename == 'boost::urls::url':
            return UrlPrinter(val)
        elif typename == 'boost::urls::url_base':
            return UrlBasePrinter(val)
        elif typename == 'boost::urls::url_view_base':
            return UrlViewBasePrinter(val)
        elif typename == 'boost::urls::detail::url_impl':
            return UrlImplPrinter(val)
        elif typename == 'boost::core::basic_string_view<char>':
            return StringViewPrinter(val)
        elif typename == 'boost::urls::pct_string_view':
            return PctStringViewPrinter(val)
        elif typename == 'boost::urls::detail::parts_base::from':
            return EnumPrinter(val)
        elif typename == 'boost::urls::scheme':
            return EnumPrinter(val)
        elif typename == 'boost::urls::host_type':
            return UrlHostTypePrinter(val)
        elif typename.startswith('boost::system::result<'):
            return SystemResultPrinter(val)
        # elif typename.startswith('boost::variant2::variant<'):
        #     return UrlHostTypePrinter(val)
        return None


def register_boost_url_printers(obj=None, dev: bool = None):
    if obj is None:
        obj = gdb
    if dev is not None:
        UrlImplPrinter.show_offset_parts = dev
    obj.pretty_printers.append(lookup_function)


if __name__ == "__main__":
    assert utils.pct_decode("Hello%20World%21") == 'Hello World!'
