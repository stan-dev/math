#!/usr/bin/env python

#
# Copyright (c) 2024 Dmitry Arkhipov (grisumbras@yandex.ru)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

import argparse
import importlib
import io
import jinja2
import json
import os.path
import re
import sys
import xml.etree.ElementTree as ET


class Nullcontext():
    def __init__(self):
        pass

    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_value, traceback):
        pass


class Access:
    public = 'public'
    protected = 'protected'
    private = 'private'

class FunctionKind:
    static = 'static'
    nonstatic = 'nonstatic'
    friend = 'friend'
    free = 'free'

class VirtualKind:
    nonvirtual = 'non-virtual'
    purevirtual = 'pure-virtual'
    virtual = 'virtual'


class Linebreak():
    pass

class PhraseContainer:
    def __init__(self, parts):
        self._parts = parts

    def __getitem__(self, pos):
        return self._parts[pos]

    def __setitem__(self, pos, val):
        self._parts[pos] = val

    def __len__(self):
        return len(self._parts)

class Phrase(PhraseContainer):
    @property
    def text(self):
        return ''.join((
            part if isinstance(part, str)
                else '' if isinstance(part, Linebreak)
                else part.text
            for part in self
        ))

class Emphasised(Phrase):
    pass

class Monospaced(Phrase):
    pass

class Strong(Phrase):
    pass

class EntityRef(Phrase):
    def __init__(self, entity, parts):
        super().__init__(parts)
        self.entity = entity

    @property
    def text(self):
        result = super().text
        if not result:
            result = self.entity.fully_qualified_name
        return result

class UrlLink(Phrase):
    def __init__(self, url, parts):
        super().__init__(parts)
        self.url = url
        assert url

    @property
    def text(self):
        result = super().text
        if not result:
            result = self.url
        return result


class Block:
    pass

class Paragraph(PhraseContainer, Block):
    def __init__(self, parts):
        super().__init__(parts)

    def __getitem__(self, pos):
        return self._parts[pos]

    def __len__(self):
        return len(self._parts)

    @property
    def text(self):
        return ''.join([
            (p if isinstance(p, str) else p.text) for p in self._parts
        ])

class List(Block):
    Arabic = '1'
    LowerLatin = 'a'
    UpperLatin = 'A'
    LowerRoman = 'i'
    UpperRoman = 'I'

    def __init__(self, kind, items):
        self.kind = kind
        self.is_ordered = kind is not None
        self._items = items

    def __getitem__(self, pos):
        return self._items[pos]

    def __len__(self):
        return len(self._items)

class ListItem(Block):
    def __init__(self, blocks):
        assert blocks
        self._blocks = blocks

    def __getitem__(self, pos):
        return self._blocks[pos]

    def __len__(self):
        return len(self.blocks)

class Section(Block):
    See = 'see'
    Returns = 'return'
    Author = 'author'
    Authors = 'authors'
    Version = 'version'
    Since = 'since'
    Date = 'date'
    Note = 'note'
    Warning = 'warning'
    Preconditions = 'pre'
    Postconditions = 'post'
    Copyright = 'copyright'
    Invariants = 'invariant'
    Remarks = 'remark'
    Attention = 'attention'
    Custom = 'par'
    RCS = 'rcs'

    def __init__(self, kind, title, blocks):
        self.kind = kind
        self.title = title
        self._blocks = blocks

    def __getitem__(self, pos):
        return self._blocks[pos]

    def __len__(self):
        return len(self._blocks)

class ParameterList(Block):
    Parameters = 'param'
    ReturnValues = 'retval'
    Exceptions = 'exception'
    TemplateParameters = 'templateparam'

    def __init__(self, kind, items):
        self.kind = kind
        self._items = items

    def __getitem__(self, pos):
        return self._items[pos]

    def __len__(self):
        return len(self._items)

class ParameterDescription(Block):
    def __init__(self, description, params):
        self.description = description
        self._params = params

    def __getitem__(self, pos):
        return self._params[pos]

    def __len__(self):
        return len(self._params)

class ParameterItem(Block):
    def __init__(self, type, name, direction):
        self.type = type
        self.name = name
        self.direction = direction

    @property
    def is_in(self):
        return self.direction in ('in', 'inout')

    @property
    def is_out(self):
        return self.direction in ('out', 'inout')

class CodeBlock(Block):
    def __init__(self, lines):
        self._lines = lines

    def __getitem__(self, pos):
        return self._lines[pos]

    def __len__(self):
        return len(self._lines)

class Table(Block):
    def __init__(self, cols, rows, caption=None, width=None):
        self.cols = cols
        self.width = width
        self.caption = caption
        self._rows = rows

    def __getitem__(self, pos):
        return self._rows[pos]

    def __len__(self):
        return len(self._rows)

class Cell(Block):
    def __init__(
        self, blocks,
        col_span=1, row_span=1, is_header=False, horizontal_align=None,
        vertical_align=None, width=None, role=None
    ):
        self._blocks = blocks

        self.col_span = int(col_span or 1)
        self.row_span = int(row_span or 1)
        self.is_header = is_header or False
        self.horizontal_align = horizontal_align
        self.vertical_align = vertical_align
        self.width = width
        self.role = role

    def __getitem__(self, pos):
        return self._blocks[pos]

    def __len__(self):
        return len(self._blocks)

def make_blocks(element, index):
    if element is None:
        return []

    result = []
    cur_para = []

    def finish_paragraph(cur_para):
        if cur_para:
            if isinstance(cur_para[0], str):
                cur_para[0] = cur_para[0].lstrip()
                if not cur_para[0]:
                    cur_para = cur_para[1:]
        if cur_para:
            if isinstance(cur_para[-1], str):
                cur_para[-1] = cur_para[-1].rstrip()
                if not cur_para[-1]:
                    cur_para = cur_para[:-1]

        # spaces after linebreaks usually cause issues
        for n in range(1, len(cur_para)):
            if (isinstance(cur_para[n - 1], Linebreak)
                    and isinstance(cur_para[n], str)):
                cur_para[n] = cur_para[n].lstrip()

        if cur_para:
            result.append(Paragraph(cur_para))

        return []

    if element.text:
        cur_para.append(remove_endlines(element.text))

    for child in element:
        func = {
            'itemizedlist': make_list,
            'simplesect': make_section,
            'programlisting': make_codeblock,
            'parameterlist': make_parameters,
            'table': make_table,
        }.get(child.tag)
        if func:
            cur_para = finish_paragraph(cur_para)
            result.append(func(child, index))
        elif child.tag == 'para':
            cur_para = finish_paragraph(cur_para)
            result.extend(make_blocks(child, index))
        else:
            cur_para.append(make_phrase(child, index))
        if child.tail:
            cur_para.append(remove_endlines(child.tail))

    finish_paragraph(cur_para)
    return result

def make_list(element, index):
    items = []
    for child in element:
        assert child.tag == 'listitem'
        items.append(make_blocks(child, index))
    return List(element.get('type'), items)

def make_parameters(element, index):
    result = []
    descr = None
    kind = element.get('kind')
    for descr_block in element:
        assert descr_block.tag == 'parameteritem'
        descr = None
        params = []
        for item in descr_block:
            if item.tag == 'parameterdescription':
                assert descr == None
                descr = item
                continue

            assert item.tag == 'parameternamelist'

            name = item.find('parametername')
            direction = name.get('direction') if name is not None else None
            params.append(
                ParameterItem(
                    text_with_refs(item.find('parametertype'), index),
                    text_with_refs(name, index),
                    direction))

        assert descr is not None
        result.append(ParameterDescription(make_blocks(descr, index), params))
    return ParameterList(kind, result)

def make_section(element, index):
    title = None
    if len(element) and element[0].tag == 'title':
        title = phrase_content(element[0], index)
    title = Paragraph(title or [])

    kind = element.get('kind')

    parts = []
    for child in element:
        if child.tag == 'para':
            parts.extend(make_blocks(child, index))

    return Section(kind, title, parts)

def make_codeblock(element, index):
    lines = []
    for line in element:
        assert line.tag == 'codeline'
        text = ''
        for hl in line:
            assert hl.tag == 'highlight'
            if hl.text:
                text += hl.text

            for part in hl:
                if part.tag == 'sp':
                    text += ' '
                elif part.tag == 'ref':
                    text += part.text or ''
                if part.tail:
                    text += part.tail
            if hl.tail:
                text += hl.tail
        lines.append(text)

    return CodeBlock(lines)

def make_table(element, index):
    cols = element.get('cols')

    caption = None
    if len(element) and element[0].tag == 'caption':
        caption = phrase_content(element[0], index)
        caption = Paragraph(caption or [])

    rows = []
    for row in element[(1 if caption else 0):]:
        assert row.tag == 'row'
        cells = []
        for cell in row:
            cells.append(Cell(
                make_blocks(cell, index),
                col_span=cell.get('colspan'),
                row_span=cell.get('rowspan'),
                is_header=cell.get('thead'),
                horizontal_align=cell.get('align'),
                vertical_align=cell.get('valign'),
                width=cell.get('width'),
                role=cell.get('class'),
            ))
        rows.append(cells)

    return Table(cols, rows, caption)

def phrase_content(element, index, allow_missing_refs=False):
    if element is None:
        return []

    result = []
    if element.text:
        result.append(remove_endlines(element.text))

    for child in element:
        result.append(
            make_phrase(child, index, allow_missing_refs=allow_missing_refs))
        if child.tail:
            result.append(remove_endlines(child.tail))

    return result

def make_phrase(element, index, allow_missing_refs=False):
    func = {
        'bold': make_strong,
        'computeroutput': make_monospaced,
        'verbatim': make_monospaced,
        'emphasis': make_emphasised,
        'ulink': make_url_link,
        'linebreak': make_linebreak,
        'ref': make_entity_reference,
    }[element.tag]
    return func(element, index, allow_missing_refs=allow_missing_refs)

def make_strong(element, index, allow_missing_refs=False):
    return Strong(
        phrase_content(element, index, allow_missing_refs=allow_missing_refs))

def make_monospaced(element, index, allow_missing_refs=False):
    return Monospaced(
        phrase_content(element, index, allow_missing_refs=allow_missing_refs))

def make_emphasised(element, index, allow_missing_refs=False):
    return Emphasised(
        phrase_content(element, index, allow_missing_refs=allow_missing_refs))

def make_url_link(element, index, allow_missing_refs=False):
    return UrlLink(
        element.get('url'),
        phrase_content(element, index, allow_missing_refs=allow_missing_refs))

def make_linebreak(element, index, allow_missing_refs=None):
    return Linebreak()

def make_entity_reference(
    element, index, refid=None, allow_missing_refs=False
):
    refid = refid or element.get('refid')
    assert refid

    target = index.get(refid)
    if target:
        return EntityRef(
            target,
            phrase_content(
                element, index, allow_missing_refs=allow_missing_refs))

    if allow_missing_refs:
        return Phrase(phrase_content(element, index, allow_missing_refs=True))

def text_with_refs(element, index):
    return Phrase(phrase_content(element, index, allow_missing_refs=True))

def resolve_type(element, index):
    result = text_with_refs(element, index)
    if (
        result
        and isinstance(result[0], str)
        and result[0].startswith('constexpr')
    ):
        result[0] = result[0][len('constexpr'):].lstrip()
    return result

_chartable = {
    ord('\r'): None,
    ord('\n'): None,
}
def remove_endlines(s):
    return s.translate(_chartable)

_noexcept_pattern = re.compile(r'(?<=\bnoexcept\()')
def parse_noexcept_condition(argstring):
    match = _noexcept_pattern.search(argstring)
    if match:
        parens = 1
        start = match.start(0)
        end = -1
        for i in range(start, len(argstring)):
            if not parens:
                end = i
                break
            if argstring[i] == '(':
                parens += 1
            elif argstring[i]== ')':
                parens -= 1
        assert parens == 0
        return argstring[start:end]

class Location():
    def __init__(self, elem):
        self.file = elem.get('file')

        self.line = elem.get('line')
        if self.line:
            self.line = int(self.line)

        self.column = elem.get('column')
        if self.column :
            self.column = int(self.column)


class Entity():
    access = Access.public

    def __init__(self, element, scope, index=dict()):
        self.id = element.get('id')
        assert self.id

        self.scope = scope

        self.name = ''.join( element.find(self.nametag).itertext() )

        self.groups = []

        loc = element.find('location')
        self._location = Location(loc) if (loc is not None) else None

        self._brief = element.find('briefdescription')
        self._description = element.find('detaileddescription')

        self.index = index
        index[self.id] = self

    @property
    def location(self):
        return (
            self._location
            or (self.scope.location if self.scope else None)
        )

    @property
    def fully_qualified_name(self):
        return '::'.join((c.name for c in self.path))

    @property
    def path(self):
        result = self.scope.path if self.scope else []
        result.append(self)
        return result

    def lookup(self, qname):
        name_parts = qname.split('::')
        if not name_parts:
            return

        scope = self
        while scope is not None:
            if scope.name == name_parts[0]:
                break
            else:
                found = None
                for entity in scope.members.values():
                    if entity.name == name_parts[0]:
                        found = entity
                        break
                if found is None:
                    scope = scope.scope
                else:
                    break
        if not scope:
            return

        for part in name_parts:
            if scope.name != part:
                found = None
                for entity in scope.members.values():
                    if entity.name == part:
                        found = entity
                        break
                scope = found
                if found is None:
                    break

        return scope

    def resolve_references(self):
        self.brief = make_blocks(self._brief, self.index)
        delattr(self, '_brief')

        self.description = make_blocks(self._description, self.index)
        delattr(self, '_description')

    def update_scopes(self):
        pass

    def __lt__(self, other):
        return self.name < other.name


class Compound(Entity):
    nametag = 'compoundname'

    def __init__(self, element, scope, index=dict()):
        super().__init__(element, scope, index)

        self.members = {}
        self._nested = []
        for section in element:
            if section.tag in ('innerclass', 'innernamespace', 'innergroup'):
                self._nested.append((
                    section.get('prot', Access.public),
                    section.get('refid')))

    def update_scopes(self):
        super().update_scopes()

        for access, refid in self._nested:
            entity = self.index[refid]
            entity.scope = self
            entity.access = access

            assert entity.name not in self.members
            self.members[entity.name] = entity
        delattr(self, '_nested')

    def resolve_references(self):
        super().resolve_references()
        if getattr(self, '_nested', None):
            self.update_scopes()


class Member(Entity):
    nametag = 'name'

    def __init__(self, element, scope, index=dict()):
        super().__init__(element, scope, index)
        self.access = element.get('prot') or Access.public


class Group(Compound):
    nametag = 'title'

    def __init__(self, element, index=dict()):
        super().__init__(element, None, index)

    def adopt(self, entity, access):
        entity.groups.append(self)


class Templatable(Entity):
    def __init__(self, element, scope, index=dict()):
        super().__init__(element, scope, index)
        self._template_parameters = element.find('templateparamlist')
        self.is_specialization = (
            (self.name.find('<') > 0) and (self.name.find('>') > 0))

    def resolve_references(self):
        super().resolve_references()
        params = (
            self._template_parameters
            if self._template_parameters is not None
                and len(self._template_parameters)
            else []
        )
        self.template_parameters = [Parameter(elem, self) for elem in params]
        delattr(self, '_template_parameters')


class Type(Templatable):
    declarator = None
    objects = []


class Scope(Entity):
    def __init__(self, element, scope, index=dict()):
        super().__init__(element, scope, index)

        nesting = 0
        name = ''
        colon = False
        for c in self.name:
            if colon:
                colon = False
                if c == ':':
                    name = ''
                    continue
                else:
                    name += ':'

            if nesting:
                if c == '<':
                    nesting += 1
                elif c == '>':
                    nesting -= 1
                name += c
            elif c == ':':
                colon = True
            else:
                if c == '<':
                    nesting += 1
                name += c
        if colon:
            name.append(':')
        self.name = name

        self.members = dict()
        for section in element:
            if not section.tag == 'sectiondef':
                continue

            for member_def in section:
                assert member_def.tag == 'memberdef'

                kind = member_def.get('kind')
                if (kind == 'friend'
                    and member_def.find('type').text == 'class'):
                    # ignore friend classes
                    continue

                factory = {
                    'function': OverloadSet.create,
                    'friend': OverloadSet.create,
                    'variable': Variable,
                    'typedef': TypeAlias,
                    'enum': Enum,
                }[kind]
                member = factory(member_def, section, self, index)
                if member is None:
                    continue
                if type(member) is OverloadSet:
                    key = (member.name, member.access, member.kind)
                    self.members[key] = member
                else:
                    assert member.name not in self.members
                    self.members[member.name] = member


    def resolve_references(self):
        super().resolve_references()
        bad_keys = []
        for key, member in self.members.items():
            if isinstance(member, OverloadSet):
                good_funcs = [
                    func for func in member
                    if self.index.get(func.id) is func
                ]
                if good_funcs:
                    member.funcs = good_funcs
                else:
                    bad_keys.append(key)
            elif self.index[member.id] is not member:
                bad_keys.append(key)

        for key in bad_keys:
            del self.members[key]

class Namespace(Scope, Compound):
    declarator = 'namespace'

    def __init__(self, element, index=dict()):
        super().__init__(element, None, index)

    def resolve_references(self):
        super().resolve_references()
        for member in self.members.values():
            if isinstance(member, OverloadSet):
                for func in member:
                    func.is_free = True


class Class(Scope, Compound, Type):
    declarator = 'class'

    def __init__(self, element, index=dict()):
        super().__init__(element, None, index)
        self._bases = element.findall('basecompoundref')

    def resolve_references(self):
        super().resolve_references()
        self.bases = [Generalization(entry, self) for entry in self._bases]
        delattr(self, '_bases')


class Generalization():
    def __init__(self, element, derived):
        self.is_virtual = element.get('virt') == 'virtual'
        self.access = element.get('prot')
        assert self.access

        refid = element.get('refid')
        if not refid:
            entity = derived.lookup( ''.join(element.itertext()).strip() )
            if entity is not None:
                refid = entity.id

        if refid:
            self.base = Phrase(
                [make_entity_reference(element, derived.index, refid)])
        else:
            self.base = text_with_refs(element, derived.index)


class Struct(Class):
    declarator = 'struct'


class Union(Class):
    declarator = 'union'


class Enum(Scope, Type, Member):
    def __init__(self, element, section, parent, index=dict()):
        super().__init__(element, parent, index)
        self.is_scoped = element.get('strong') == 'yes'
        self._underlying_type = element.find('type')

        self.objects = []
        for child in element.findall('enumvalue'):
            enumerator = Enumerator(child, section, self, index)
            self.objects.append(enumerator)
            assert enumerator.name not in enumerator.scope.members
            enumerator.scope.members[enumerator.name] = enumerator

    @property
    def declarator(self):
        return 'enum class' if self.is_scoped else 'enum'

    def resolve_references(self):
        super().resolve_references()
        self.underlying_type = text_with_refs(
            self._underlying_type, self.index)
        delattr(self, '_underlying_type')


class Value(Member, Templatable):
    def __init__(self, element, parent, index=dict()):
        super().__init__(element, parent, index)
        self.is_static = element.get('static') == 'yes'
        self.is_constexpr = element.get('constexpr') == 'yes'
        self.is_volatile = element.get('volatile') == 'yes'
        self.is_const = (
            element.get('const') == 'yes'
            or element.get('mutable') == 'no'
            or False)
        self.is_inline = element.get('inline') == 'yes'


class Function(Value):
    def __init__(self, element, section, parent, index=dict()):
        super().__init__(element, parent, index)
        self.is_explicit = element.get('explicit') == 'yes'
        self.refqual = element.get('refqual')
        self.virtual_kind = element.get('virt', VirtualKind.nonvirtual)
        self.is_friend = section.get('kind') == 'friend'
        self.is_free = section.get('kind') == 'related'
        self.is_constructor = self.name == parent.name
        self.is_destructor = self.name == '~' + parent.name
        noexcept = element.get('noexcept')
        if noexcept:
            self.is_noexcept = noexcept == 'yes'
        else:
            self.is_noexcept = self.is_destructor

        args = element.find('argsstring').text or ''
        self.is_deleted = args.endswith('=delete')
        self.is_defaulted = args.endswith('=default')

        self.overload_set = None

        self._return_type = element.find('type')
        assert (
            self.is_constructor
            or self.is_destructor
            or self._return_type is not None)

        self._parameters = element.findall('param')

        self.noexcept_condition = None
        if self.is_noexcept:
            self.noexcept_condition = element.get('noexceptexpression')
            if not self.noexcept_condition:
                self.noexcept_condition = parse_noexcept_condition(args)

    @property
    def kind(self):
        if self.is_friend:
            return FunctionKind.friend
        if self.is_free:
            return FunctionKind.free
        if self.is_static:
            return FunctionKind.static
        return FunctionKind.nonstatic

    @property
    def is_sole_overload(self):
        return len(self.overload_set) == 1

    @property
    def overload_index(self):
        if self.overload_set is None or self.is_sole_overload:
            return -1

        for n, overload in enumerate(self.overload_set):
            if self == overload:
                return n

        return -1

    def resolve_references(self):
        super().resolve_references()

        self.return_type = resolve_type(self._return_type, self.index)
        delattr(self, '_return_type')

        self.parameters = [Parameter(elem, self) for elem in self._parameters]
        delattr(self, '_parameters')

    def __lt__(self, other):
        if not isinstance(other, Function):
            return self.name < other.name

        if self.name == other.name:
            if self.scope == other.scope:
                return self.overload_index < other.overload_index
            return self.scope < other.scope

        # if self.is_constructor:
        #     return True
        # if other.is_constructor:
        #     return False
        # if self.is_destructor:
        #     return True
        # if other.is_destructor:
        #     return False
        return self.name < other.name


class Parameter():
    def __init__(self, element, parent):
        self.type = text_with_refs(element.find('type'), parent.index)
        self.default_value = text_with_refs(
            element.find('defval'), parent.index)

        self.description = element.find('briefdescription')
        if self.description is not None:
            self.description = make_blocks(self.description, parent)
        else:
            self.description = []

        self.name = element.find('declname')
        if self.name is not None:
            self.name = self.name.text

        self.array = text_with_refs(element.find('array'), parent.index)
        if self.array:
            assert isinstance(self.type[-1], str)
            assert self.type[-1].endswith('(&)')
            self.type[-1] = self.type[-1][:-3]

        self.args = text_with_refs(element.find('argsstring'), parent.index)
        if self.args:
            assert isinstance(self.type[-1], str)
            assert isinstance(self.args[0], str)
            assert self.type[-1].endswith('(*')
            assert self.args[0].startswith(')(')
            self.type[-1] = self.type[-1][:-2]
            self.args[0] = self.args[0][1:]


class OverloadSet():
    @staticmethod
    def create(element, section, parent, index):
        if (index.get( element.get('id') )
            and (section.get('kind') != 'related'
                 or not isinstance(parent, Class))
        ):
            return None

        func = Function(element, section, parent, index)

        key = (func.name, func.access, func.kind)
        if key in parent.members:
            overload_set = parent.members[key]
            overload_set.append(func)
        else:
            overload_set = OverloadSet([func])
        func.overload_set = overload_set
        return overload_set

    def __init__(self, funcs):
        self.funcs = funcs
        assert len(funcs)
        self._resort()

    def append(self, func):
        self.funcs.append(func)
        self._resort()

    def __getattr__(self, name):
        return getattr(self.funcs[0], name)

    @property
    def brief(self):
        return [func.brief for func in self.funcs]

    def __getitem__(self, pos):
        return self.funcs[pos]

    def __len__(self):
        return len(self.funcs)

    def __lt__(self, other):
        if isinstance(other, OverloadSet):
            return self.funcs[0] < other.funcs[0]
        return self.name < other.name

    def _resort(self):
        funcs_backwards = []
        briefs_backwards = []
        for func in self.funcs:
            brief = (
                ''.join(func._brief.itertext())
                if func._brief is not None
                else ''
            )
            try:
                n = briefs_backwards.index(brief)
            except:
                n = 0
            briefs_backwards.insert(n, brief)
            funcs_backwards.insert(n, func)
        funcs_backwards.reverse()
        self.funcs = funcs_backwards


class Variable(Value):
    def __init__(self, element, section, parent, index=dict()):
        super().__init__(element, parent, index)
        self._value = element.find('initializer')
        self._type = element.find('type')
        self._args = element.find('argsstring')

    def resolve_references(self):
        super().resolve_references()

        self.value = text_with_refs(self._value, self.index)
        delattr(self, '_value')

        self.type = resolve_type(self._type, self.index)
        delattr(self, '_type')

        self.args = text_with_refs(self._args, self.index)
        delattr(self, '_args')
        if (self.args
            and self.type[-1].endswith('(*')
            and self.args[0].startswith(')(')
        ):
            self.type[-1] = self.type[-1][:-2]
            self.args[0] = self.args[0][1:]


class Enumerator(Variable):
    def __init__(self, element, section, parent, index=dict()):
        super().__init__(element, section, parent, index)
        self.is_constexpr = True
        self.is_const = True
        self.is_static = True
        self.enum = parent
        self.scope = parent if parent.is_scoped else parent.scope
        assert self.scope

    def resolve_references(self):
        super().resolve_references()
        self.type = Phrase([EntityRef(self.enum, [self.enum.name])])


class TypeAlias(Member, Type):
    declarator = 'using'

    def __init__(self, element, section, parent, index=dict()):
        super().__init__(element, parent, index)

        self._aliased = element.find('type')
        assert self._aliased is not None

    def resolve_references(self):
        super().resolve_references()
        self.aliased = text_with_refs(self._aliased, self.index)
        delattr(self, '_aliased')


class AcceptOneorNone(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        if kwargs.get('nargs') is not None:
            raise ValueError("nargs not allowed")
        if kwargs.get('const') is not None:
            raise ValueError("const not allowed")
        if kwargs.get('default') is not None:
            raise ValueError("default not allowed")

        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if getattr(namespace, self.dest) is not None:
            raise argparse.ArgumentError(self, "multiple values")

        setattr(namespace, self.dest, values)

def parse_args(args):
    parser = argparse.ArgumentParser(
        prog=args[0],
        description='Produces API reference in QuickBook markup')
    parser.add_argument(
        '-i',
        '--input',
        action=AcceptOneorNone,
        help='Doxygen XML index file; STDIN by default')
    parser.add_argument(
        '-o',
        '--output',
        action=AcceptOneorNone,
        help='Output file; STDOUT by default')
    parser.add_argument(
        '-c', '--config',
        action='append',
        default=[],
        help='Configuration files')
    parser.add_argument(
        '-T', '--template',
        action=AcceptOneorNone,
        help='Jinja2 template to use for output')
    parser.add_argument(
        '-E', '--extension',
        action='append',
        default=[],
        help='Extension module')
    parser.add_argument(
        '-I', '--include',
        action='append',
        default=[],
        help='Directory with template partials')
    parser.add_argument(
        '-D', '--directory',
        action=AcceptOneorNone,
        help=(
            'Directory with additional data files; '
            'by default INPUT parent directory if that is provided, '
            'otherwise PWD'))
    return parser.parse_args(args[1:])

def open_input(stdin, args, cwd):
    data_dir = args.directory
    if args.input:
        file = open(args.input, 'r', encoding='utf-8')
        ctx = file
        data_dir = data_dir or os.path.dirname(args.input)
    else:
        file = stdin
        ctx = Nullcontext()
        data_dir = data_dir or cwd
    return (file, ctx, data_dir)

def open_output(stdout, args):
    if args.output:
        file = open(args.output, 'w', encoding='utf-8')
        ctx = file
    else:
        file = stdout
        ctx = Nullcontext()
    return (file, ctx)

def load_configs(args):
    result = {
        'include_private': False,
        'legacy_behavior': True,
    }
    for file_name in args.config:
        with open(file_name, 'r', encoding='utf-8') as file:
            result.update( json.load(file) )
    return result

def collect_compound_refs(file):
    tree = ET.parse(file)
    root = tree.getroot()
    for child in root:
        if child.tag != 'compound':
            continue

        kind = child.get('kind')
        assert kind
        if kind in ('file', 'dir'):
            continue

        refid = child.get('refid')
        assert refid
        yield refid

def collect_data(parent_dir, refs):
    result = dict()
    for refid in refs:
        file_name = os.path.join(parent_dir, refid) + '.xml'
        with open(file_name, 'r', encoding='utf-8') as file:
            tree = ET.parse(file)
            root = tree.getroot()
            assert len(root) == 1

            element = root[0]
            assert element.tag == 'compounddef'

            factory = {
                'class': Class,
                'namespace': Namespace,
                'struct': Struct,
                'union': Union,
                'group': Group
            }.get(element.get('kind'))
            if not factory:
                continue
            factory(element, result)

    for entity in result.values():
        assert entity is not None
        entity.update_scopes()

    for entity in result.values():
        entity.resolve_references();

    return result

def docca_include_dir(script):
    return os.path.join(os.path.dirname(script), 'include')

def template_file_name(includes, args):
    return (args.template
         or os.path.join(includes, 'docca/quickbook.jinja2'))

def collect_include_dirs(template, include_dir, args):
    result =  [os.path.dirname(template), include_dir]
    result.extend(args.include)
    return result

def construct_environment(loader, config):
    env = jinja2.Environment(
        loader=loader,
        autoescape=False,
        undefined=jinja2.StrictUndefined,
        extensions=[
            'jinja2.ext.do',
            'jinja2.ext.loopcontrols',
        ],
    )

    env.globals['Access'] = Access
    env.globals['FunctionKind'] = FunctionKind
    env.globals['VirtualKind'] = VirtualKind
    env.globals['Section'] = Section
    env.globals['ParameterList'] = ParameterList
    env.globals['Config'] = config
    env.globals['re'] = re

    env.tests['Entity'] = lambda x: isinstance(x, Entity)
    env.tests['Templatable'] = lambda x: isinstance(x, Templatable)
    env.tests['Type'] = lambda x: isinstance(x, Type)
    env.tests['Scope'] = lambda x: isinstance(x, Scope)
    env.tests['Namespace'] = lambda x: isinstance(x, Namespace)
    env.tests['TypeAlias'] = lambda x: isinstance(x, TypeAlias)
    env.tests['Class'] = lambda x: isinstance(x, Class)
    env.tests['Struct'] = lambda x: isinstance(x, Struct)
    env.tests['Union'] = lambda x: isinstance(x, Union)
    env.tests['Enum'] = lambda x: isinstance(x, Enum)
    env.tests['Value'] = lambda x: isinstance(x, Value)
    env.tests['Variable'] = lambda x: isinstance(x, Variable)
    env.tests['Enumerator'] = lambda x: isinstance(x, Enumerator)
    env.tests['Function'] = lambda x: isinstance(x, Function)
    env.tests['OverloadSet'] = lambda x: isinstance(x, OverloadSet)
    env.tests['Parameter'] = lambda x: isinstance(x, Parameter)

    env.tests['Phrase'] = lambda x: isinstance(x, Phrase)
    env.tests['Linebreak'] = lambda x: isinstance(x, Linebreak)
    env.tests['Emphasised'] = lambda x: isinstance(x, Emphasised)
    env.tests['Strong'] = lambda x: isinstance(x, Strong)
    env.tests['Monospaced'] = lambda x: isinstance(x, Monospaced)
    env.tests['EntityRef'] = lambda x: isinstance(x, EntityRef)
    env.tests['UrlLink'] = lambda x: isinstance(x, UrlLink)

    env.tests['Block'] = lambda x: isinstance(x, Block)
    env.tests['Paragraph'] = lambda x: isinstance(x, Paragraph)
    env.tests['List'] = lambda x: isinstance(x, List)
    env.tests['ListItem'] = lambda x: isinstance(x, ListItem)
    env.tests['Section'] = lambda x: isinstance(x, Section)
    env.tests['CodeBlock'] = lambda x: isinstance(x, CodeBlock)
    env.tests['Table'] = lambda x: isinstance(x, Table)
    env.tests['Cell'] = lambda x: isinstance(x, Cell)
    env.tests['ParameterList'] = lambda x: isinstance(x, ParameterList)
    env.tests['ParameterDescription'] = lambda x: isinstance(x, ParameterDescription)
    env.tests['ParameterItem'] = lambda x: isinstance(x, ParameterItem)

    return env

def load_extensions(files):
    result = []
    counter = 0
    for file in files:
        module = None
        if not os.path.exists(file):
            raise RuntimeError('Could not find module %s' % file)

        name = 'docca._ext' + str(counter)
        spec = importlib.util.spec_from_file_location(name, file)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        result.append(module)
        counter += 1

    return result

def install_extensions(env, exts):
    for ext in exts:
        ext.install_docca_extension(env)
    return env

def render(env, file_name, output, data):
    template = env.get_template(os.path.basename(file_name))
    template.stream(entities=data).dump(output)

def main(args, stdin, stdout, script):
    args = parse_args(args)

    file, ctx, data_dir = open_input(stdin, args, os.getcwd())
    with ctx:
        refs = list(collect_compound_refs(file))
    data = collect_data(data_dir, refs)

    config = load_configs(args)

    file, ctx = open_output(stdout, args)
    with ctx:
        include_dir = docca_include_dir(script)
        template = template_file_name(include_dir, args)

        include_dirs = collect_include_dirs(template, include_dir, args)

        env = construct_environment(
            jinja2.FileSystemLoader(include_dirs), config)

        exts = load_extensions(args.extension)
        env = install_extensions(env, exts)

        render(env, template, file, data)

if __name__ == '__main__':
    main(sys.argv, sys.stdin, sys.stdout, os.path.realpath(__file__))
