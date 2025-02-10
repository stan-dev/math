#
# Copyright (c) 2024 Dmitry Arkhipov (grisumbras@yandex.ru)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

import argparse
import contextlib
import docca
import io
import jinja2
import jinja2.ext
import os
import pytest
import re
import types

from docca_test_helpers import (
    MockXmlElem,
    make_elem,
)

def collect_paragraphs(pars):
    return ''.join([p.text for p in pars])

def test_blocks():
    assert not docca.make_blocks(make_elem({}), None)

    blocks = docca.make_blocks(
        make_elem({'items': [
            'some text',
            { 'tag': 'para', 'items': ['more text'] },
            { 'tag': 'bold', 'items': ['bold text'] },
            ' final text',
        ]}),
        None)
    assert len(blocks) == 3

    assert isinstance(blocks[0], docca.Paragraph)
    assert len(blocks[0]) == 1
    assert blocks[0].text == 'some text'
    assert blocks[0][0] == 'some text'

    assert isinstance(blocks[1], docca.Paragraph)
    assert len(blocks[1]) == 1
    assert blocks[1].text == 'more text'
    assert blocks[1][0] == 'more text'

    assert isinstance(blocks[2], docca.Paragraph)
    assert len(blocks[2]) == 2
    assert blocks[2].text == 'bold text final text'
    assert blocks[2][0].text == 'bold text'
    assert blocks[2][1] == ' final text'

    blocks = docca.make_blocks(
        make_elem({'items': [
            { 'tag': 'simplesect'},
        ]}),
        None)
    assert len(blocks) == 1
    assert isinstance(blocks[0], docca.Section)
    assert isinstance(blocks[0].title, docca.Paragraph)
    assert not blocks[0].title
    assert not blocks[0]
    assert not blocks[0].kind

    blocks = docca.make_blocks(
        make_elem({'items': [
            { 'tag': 'simplesect', 'kind': 'some kind' },
        ]}),
        None)
    assert len(blocks) == 1
    assert isinstance(blocks[0], docca.Section)
    assert isinstance(blocks[0].title, docca.Paragraph)
    assert not blocks[0].title
    assert not blocks[0]
    assert blocks[0].kind == 'some kind'

    blocks = docca.make_blocks(
        make_elem({'items': [
            {
                'tag': 'simplesect',
                'items': [
                    {
                        'tag': 'title',
                        'items': [
                            'some text ',
                            { 'tag': 'bold', 'items': ['bold text'] },
                        ],
                    },
                ]
            },
        ]}),
        None)
    assert not blocks[0]
    assert blocks[0].title
    assert blocks[0].title.text == 'some text bold text'
    assert blocks[0].title[0] == 'some text '
    assert blocks[0].title[1].text == 'bold text'

    blocks = docca.make_blocks(
        make_elem({'items': [
            {
                'tag': 'simplesect',
                'items': [
                    {
                        'tag': 'title',
                        'items': [
                            'some text ',
                            { 'tag': 'bold', 'items': ['bold text'] },
                        ],
                    },
                    { 'tag': 'para', 'items': ['paragraph one'] },
                    { 'tag': 'para', 'items': ['paragraph two'] },
                ]
            },
        ]}),
        None)
    assert len(blocks[0]) == 2
    assert isinstance(blocks[0][0], docca.Paragraph)
    assert blocks[0][0].text == 'paragraph one'
    assert isinstance(blocks[0][1], docca.Paragraph)
    assert blocks[0][1].text == 'paragraph two'

    blocks = docca.make_blocks(
        make_elem({'items': [
            {
                'tag': 'programlisting',
                'items': [
                    {
                        'tag': 'codeline',
                        'items': [
                            { 'tag': 'highlight', 'items': ['int n = 0;'] },
                        ],
                    },
                    {
                        'tag': 'codeline',
                        'items': [
                            {
                                'tag': 'highlight',
                                'items': [
                                    'int',
                                    { 'tag': 'sp' },
                                    { 'tag': 'ref', 'items': ['m'] },
                                    ' = 0;',
                                ],
                            },
                        ],
                    },
                ]
            },
        ]}),
        None)
    assert len(blocks) == 1
    assert isinstance(blocks[0], docca.CodeBlock)
    assert len(blocks[0]) == 2
    assert blocks[0][0] == 'int n = 0;'
    assert blocks[0][1] == 'int m = 0;'

    blocks = docca.make_blocks(
        make_elem({'items': [
            { 'tag': 'parameterlist' },
        ]}),
        None)
    assert len(blocks) == 1
    assert isinstance(blocks[0], docca.ParameterList)
    assert not blocks[0]
    assert not blocks[0].kind

    blocks = docca.make_blocks(
        make_elem({'items': [
            { 'tag': 'parameterlist', 'kind': 'somekind' },
        ]}),
        None)
    assert not blocks[0]
    assert blocks[0].kind == 'somekind'

    blocks = docca.make_blocks(
        make_elem({'items': [
            {
                'tag': 'parameterlist',
                'items': [
                    {
                        'tag': 'parameteritem',
                        'items': [
                            {
                                'tag': 'parameterdescription',
                                'items': ['one param kind']
                            },
                            {
                                'tag': 'parameternamelist',
                                'items': [
                                    {
                                        'tag': 'parametername',
                                        'direction': 'in',
                                        'items': ['A'],
                                    },
                                    {
                                        'tag': 'parametertype',
                                        'items': ['T'],
                                    },
                                ],
                            },
                            {
                                'tag': 'parameternamelist',
                                'items': [
                                    {
                                        'tag': 'parametername',
                                        'direction': 'out',
                                        'items': ['B'],
                                    },
                                    {
                                        'tag': 'parametertype',
                                        'items': ['U'],
                                    },
                                ],
                            },
                        ],
                    },
                ],
            },
        ]}),
        None)
    assert len(blocks[0]) == 1
    assert isinstance(blocks[0][0], docca.ParameterDescription)
    assert len(blocks[0][0].description) == 1
    assert isinstance(blocks[0][0].description[0], docca.Paragraph)
    assert blocks[0][0].description[0].text == 'one param kind'
    assert len(blocks[0][0]) == 2
    assert isinstance(blocks[0][0][0], docca.ParameterItem)
    assert isinstance(blocks[0][0][1], docca.ParameterItem)
    assert blocks[0][0][0].name.text == 'A'
    assert blocks[0][0][1].name.text == 'B'
    assert blocks[0][0][0].type.text == 'T'
    assert blocks[0][0][1].type.text == 'U'
    assert blocks[0][0][0].direction == 'in'
    assert blocks[0][0][1].direction == 'out'
    assert blocks[0][0][0].is_in
    assert not blocks[0][0][1].is_in
    assert not blocks[0][0][0].is_out
    assert blocks[0][0][1].is_out

    blocks = docca.make_blocks(
        make_elem({'items': [
            { 'tag': 'itemizedlist' },
        ]}),
        None)
    assert len(blocks) == 1
    assert isinstance(blocks[0], docca.List)
    assert not blocks[0]
    assert not blocks[0].kind
    assert not blocks[0].is_ordered

    blocks = docca.make_blocks(
        make_elem({'items': [
            { 'tag': 'itemizedlist', 'type': 'i' },
        ]}),
        None)
    assert blocks[0].kind == docca.List.LowerRoman
    assert blocks[0].is_ordered

    blocks = docca.make_blocks(
        make_elem({'items': [
            {
                'tag': 'itemizedlist',
                'items': [
                    { 'tag': 'listitem', 'items': ['item 1'] },
                    { 'tag': 'listitem', 'items': ['item 2'] },
                ],
            },
        ]}),
        None)
    assert len(blocks[0]) == 2
    assert len(blocks[0][0]) == 1
    assert isinstance(blocks[0][0][0], docca.Paragraph)
    assert blocks[0][0][0].text == 'item 1'
    assert len(blocks[0][1]) == 1
    assert isinstance(blocks[0][1][0], docca.Paragraph)
    assert blocks[0][1][0].text == 'item 2'

    blocks = docca.make_blocks(
        make_elem({'items': [
            { 'tag': 'table' },
        ]}),
        None)
    assert len(blocks) == 1
    assert isinstance(blocks[0], docca.Table)
    assert not blocks[0]
    assert not blocks[0].caption

    blocks = docca.make_blocks(
        make_elem({'items': [
            {
                'tag': 'table',
                'items': [ { 'tag': 'caption', 'items': ['a caption'] } ],
            }
        ]}),
        None)
    assert isinstance(blocks[0].caption, docca.Paragraph)
    assert blocks[0].caption.text == 'a caption'

    blocks = docca.make_blocks(
        make_elem({'items': [
            {
                'tag': 'table',
                'items': [
                    {
                        'tag': 'row',
                        'items': [
                            {
                                'tag': 'entry',
                                'items': [ { 'tag': 'para', 'items': ['1'] } ],
                            },
                            {
                                'tag': 'entry',
                                'thead': 'yes',
                                'align': 'left',
                                'valign': 'middle',
                                'width': '12',
                                'class': 'red',
                                'items': [ { 'tag': 'para', 'items': ['2'] } ],
                            },
                        ],
                    },
                    {
                        'tag': 'row',
                        'items': [
                            {
                                'tag': 'entry',
                                'colspan': '2',
                                'rowspan': '3',
                                'items': [ { 'tag': 'para', 'items': ['3'] } ],
                            },
                        ],
                    },
                ],
            }
        ]}),
        None)
    assert len(blocks[0]) == 2
    assert len(blocks[0][0]) == 2

    assert isinstance(blocks[0][0][0], docca.Cell)
    assert len(blocks[0][0][0]) == 1
    assert blocks[0][0][0][0].text == '1'
    assert blocks[0][0][0].col_span == 1
    assert blocks[0][0][0].row_span == 1
    assert not blocks[0][0][0].is_header
    assert not blocks[0][0][0].horizontal_align
    assert not blocks[0][0][0].vertical_align
    assert not blocks[0][0][0].width
    assert not blocks[0][0][0].role

    assert isinstance(blocks[0][0][1], docca.Cell)
    assert blocks[0][0][1][0].text == '2'
    assert blocks[0][0][1].is_header
    assert blocks[0][0][1].horizontal_align == 'left'
    assert blocks[0][0][1].vertical_align == 'middle'
    assert blocks[0][0][1].width == '12'
    assert blocks[0][0][1].role == 'red'

    assert len(blocks[0][1]) == 1
    assert blocks[0][1][0][0].text == '3'
    assert blocks[0][1][0].col_span == 2
    assert blocks[0][1][0].row_span == 3

def test_phrases():
    p = docca.make_phrase(
        make_elem({
            'tag': 'bold',
        }),
        None)
    assert isinstance(p, docca.Strong)
    assert p.text == ''
    assert not p

    p = docca.make_phrase(
        make_elem({
            'tag': 'bold',
            'items': ['bold text'],
        }),
        None)
    assert p.text == 'bold text'
    assert len(p) == 1
    assert p[0] == 'bold text'

    p = docca.make_phrase(
        make_elem({
            'tag': 'emphasis',
            'items': [
                'emphasised text ',
                {
                    'tag': 'bold',
                    'items': ['bold text'],
                },
            ],
        }),
        None)
    assert isinstance(p, docca.Emphasised)
    assert p.text == 'emphasised text bold text'
    assert len(p) == 2
    assert p[0] == 'emphasised text '
    assert p[1][0] == 'bold text'

    p = docca.make_phrase(
        make_elem({
            'tag': 'computeroutput',
            'items': [
                {
                    'tag': 'emphasis',
                    'items': ['emphasised text'],
                },
                ' computer output',
            ],
        }),
        None)
    assert isinstance(p, docca.Monospaced)
    assert p.text == 'emphasised text computer output'
    assert len(p) == 2
    assert p[0][0] == 'emphasised text'
    assert p[1] == ' computer output'

    p = docca.make_phrase(
        make_elem({
            'tag': 'verbatim',
            'items': [
                {
                    'tag': 'computeroutput',
                    'items': ['computer output'],
                },
                ' verbatim text',
            ],
        }),
        None)
    assert isinstance(p, docca.Monospaced)
    assert p.text == 'computer output verbatim text'
    assert len(p) == 2
    assert p[0][0] == 'computer output'
    assert p[1] == ' verbatim text'

    p = docca.make_phrase(
        make_elem({
            'tag': 'verbatim',
            'items': ['text', { 'tag': 'linebreak' } ],
        }),
        None)
    assert p.text == 'text'
    assert p[0] == 'text'
    assert isinstance(p[1], docca.Linebreak)

    p = docca.make_phrase(
        make_elem({
            'tag': 'ulink',
            'url': 'https:/example.tld/some/url',
            'items': [
                {
                    'tag': 'verbatim',
                    'items': ['text', { 'tag': 'linebreak' } ],
                },
            ],
        }),
        None)
    assert isinstance(p, docca.UrlLink)
    assert p.text == 'text'
    assert len(p) == 1
    assert len(p[0]) == 2
    assert p[0][0] == 'text'
    assert isinstance(p[0][1], docca.Linebreak)

    p = docca.make_phrase(
        make_elem({
            'tag': 'ref',
            'refid': 'foobar',
            'items': [ 'alternative text' ]
        }),
        dict(),
        allow_missing_refs=True)
    assert isinstance(p, docca.Phrase)
    assert p.text == 'alternative text'

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'foobar',
            'items': [ { 'tag': 'compoundname', 'items': ['ns'] } ],
        }))
    p = docca.make_phrase(
        make_elem({
            'tag': 'ref',
            'refid': 'foobar',
            'items': [ 'alternative text' ]
        }),
        dict(foobar=ns))
    assert isinstance(p, docca.Phrase)
    assert p.text == 'alternative text'
    assert p.entity == ns

    p = docca.make_phrase(
        make_elem({
            'tag': 'bold',
            'items': [ { 'tag': 'ref', 'refid': 'foobar' } ],
        }),
        dict(foobar=ns))
    assert isinstance(p[0], docca.EntityRef)
    assert p[0].entity == ns

    p = docca.make_phrase(
        make_elem({
            'tag': 'bold',
            'items': [ { 'tag': 'ref', 'refid': 'foobar' } ],
        }),
        dict(),
        allow_missing_refs=True)
    assert isinstance(p[0], docca.Phrase)
    assert p[0].text == ''

def test_namespace():
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyNs'] },
                { 'tag': 'briefdescription', 'items': ['a namespace'] },
                { 'tag': 'detaileddescription', 'items': ['this is a namespace'] },
                {
                    'tag': 'location',
                    'file': 'some/file.cpp',
                    'line': '100',
                    'column': '12'
                },
            ]
        }))
    ns.resolve_references()
    assert ns.declarator == 'namespace'
    assert not ns.scope
    assert ns.id == 'someid'
    assert ns.name == 'MyNs'
    assert ns.fully_qualified_name == 'MyNs'
    assert not ns.members
    assert collect_paragraphs(ns.brief) == 'a namespace'
    assert collect_paragraphs(ns.description) == 'this is a namespace'
    assert ns.location.file == 'some/file.cpp'
    assert ns.location.line == 100
    assert ns.location.column == 12

    index = dict()
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns1',
            'items': [
                { 'tag': 'compoundname', 'items': ['OtherNs::MyNs'] },
            ]
        }),
        index)
    ns.resolve_references()
    assert ns.name == 'MyNs'
    assert ns.fully_qualified_name == 'MyNs'

    ns2 = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns2',
            'items': [
                { 'tag': 'compoundname', 'items': ['OtherNs'] },
                { 'tag': 'innernamespace', 'refid': 'ns1' },
                {
                    'tag': 'location',
                    'file': 'some/file.cpp',
                    'line': '100',
                    'column': '12'
                },
            ]
        }),
        index)
    ns2.resolve_references()
    assert ns2.name == 'OtherNs'
    assert ns2.fully_qualified_name == 'OtherNs'
    assert ns2.members[ns.name] == ns
    assert ns.scope == ns2
    assert ns.fully_qualified_name == 'OtherNs::MyNs'
    assert ns.location.file == 'some/file.cpp'
    assert ns.location.line == 100
    assert ns.location.column == 12
    assert ns.lookup('MyNs') == ns
    assert ns.lookup('OtherNs') == ns2
    assert ns.lookup('OtherNs::MyNs') == ns
    assert ns2.lookup('OtherNs') == ns2
    assert ns2.lookup('MyNs') == ns
    assert ns2.lookup('OtherNs::MyNs') == ns

    index = dict()
    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'classid',
            'items': [
                { 'tag': 'compoundname', 'items': ['OtherNs::MyClass'] },
            ]
        }),
        index)
    c.resolve_references()

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns1',
            'items': [
                { 'tag': 'compoundname', 'items': ['OtherNs'] },
                { 'tag': 'innerclass', 'refid': 'classid' },
            ]
        }),
        index)
    ns.resolve_references()
    assert ns.name == 'OtherNs'
    assert ns.fully_qualified_name == 'OtherNs'
    assert c.name == 'MyClass'
    assert c.scope == ns
    assert c.fully_qualified_name == 'OtherNs::MyClass'
    assert c.lookup('MyClass') == c
    assert c.lookup('OtherNs') == ns
    assert c.lookup('OtherNs::MyClass') == c

    ns2.name = 'TopNs'
    ns2.members[ns.name] = ns
    ns.scope = ns2
    assert c.lookup('TopNs') == ns2
    assert c.lookup('TopNs::OtherNs') == ns
    assert c.lookup('TopNs::OtherNs::MyClass') == c
    assert ns2.lookup('TopNs') == ns2
    assert ns2.lookup('TopNs::OtherNs') == ns
    assert ns2.lookup('TopNs::OtherNs::MyClass') == c
    assert ns2.lookup('OtherNs') == ns
    assert ns2.lookup('OtherNs::MyClass') == c

def test_class():
    for Kind in (docca.Union, docca.Struct, docca.Class):
        c = Kind(
            make_elem({
                'tag': 'compound',
                'id': 'someid',
                'items': [
                    { 'tag': 'compoundname', 'items': ['MyClass'] },
                    { 'tag': 'briefdescription', 'items': ['a class'] },
                    {
                        'tag': 'detaileddescription',
                        'items': ['this is a class']
                    },
                ]
            }))
        c.resolve_references()
        assert c.id == 'someid'
        assert c.name == 'MyClass'
        assert not c.scope
        assert not c.members
        assert not c.bases
        assert not c.objects
        assert not c.template_parameters
        assert not c.is_specialization
        assert collect_paragraphs(c.brief) == 'a class'
        assert collect_paragraphs(c.description) == 'this is a class'

        c = Kind(
            make_elem({
                'tag': 'compound',
                'id': 'someid',
                'items': [
                    { 'tag': 'compoundname', 'items': ['my_ns::MyClass'] },
                    { 'tag': 'briefdescription', 'items': ['a class'] },
                    {
                        'tag': 'detaileddescription',
                        'items': ['this is a class']
                    },
                ]
            }))
        c.resolve_references()
        assert c.name == 'MyClass'

    c = docca.Struct(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [ { 'tag': 'compoundname', 'items': ['MyClass'] } ]
        }))
    c.resolve_references()
    assert c.declarator == 'struct'

    c = docca.Union(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [ { 'tag': 'compoundname', 'items': ['MyClass'] } ]
        }))
    c.resolve_references()
    assert c.declarator == 'union'

    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [ { 'tag': 'compoundname', 'items': ['MyClass'] } ]
        }))
    c.resolve_references()
    assert c.declarator == 'class'

    index = dict(someid=c)
    d = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'otherid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyClass'] },
                {
                    'tag': 'basecompoundref',
                    'virt': 'non-virtual',
                    'prot': 'public',
                    'refid': 'someid',
                },
                {
                    'tag': 'basecompoundref',
                    'virt': 'virtual',
                    'prot': 'private',
                    'refid': 'someid',
                },
                {
                    'tag': 'basecompoundref',
                    'virt': 'virtual',
                    'prot': 'protected',
                    'items': ['external::type'],
                },
            ]
        }),
        index)
    d.resolve_references()
    assert len(d.bases) == 3

    assert d.bases[0].access == 'public'
    assert not d.bases[0].is_virtual
    assert isinstance(d.bases[0].base, docca.Phrase)
    assert isinstance(d.bases[0].base[0], docca.EntityRef)
    assert d.bases[0].base[0].entity == c

    assert d.bases[1].access == 'private'
    assert d.bases[1].is_virtual
    assert d.bases[1].base[0].entity == c

    assert d.bases[2].access == 'protected'
    assert d.bases[2].is_virtual
    assert d.bases[2].base.text == 'external::type'

def test_enum():
    e = docca.Enum(
        make_elem({
            'tag': 'compound',
            'id': 'enum1',
            'items': [
                { 'tag': 'name', 'items': ['MyEnum'] },
                { 'tag': 'briefdescription', 'items': ['an enum'] },
                { 'tag': 'detaileddescription', 'items': ['this is an enum'] },
            ]
        }),
        None,
        None)
    e.resolve_references()
    assert e.declarator == 'enum'
    assert e.id == 'enum1'
    assert e.name == 'MyEnum'
    assert e.access == 'public'
    assert collect_paragraphs(e.brief) == 'an enum'
    assert collect_paragraphs(e.description) == 'this is an enum'
    assert not e.scope
    assert not e.is_scoped
    assert not e.objects
    assert not e.underlying_type.text
    assert not e.members

    e = docca.Enum(
        make_elem({
            'tag': 'compound',
            'id': 'enum1',
            'strong': 'yes',
            'items': [ { 'tag': 'name', 'items': ['MyEnum'] } ]
        }),
        None,
        None)
    e.resolve_references()
    assert e.is_scoped == True
    assert e.declarator == 'enum class'

    e = docca.Enum(
        make_elem({
            'tag': 'compound',
            'id': 'enum1',
            'strong': 'yes',
            'items': [
                { 'tag': 'name', 'items': ['MyEnum'] },
                { 'tag': 'type', 'items': ['std::size_t'] },
            ]
        }),
        None,
        None)
    e.resolve_references()
    assert e.underlying_type.text == 'std::size_t'

    e = docca.Enum(
        make_elem({
            'tag': 'compound',
            'id': 'enum1',
            'strong': 'yes',
            'items': [
                { 'tag': 'name', 'items': ['MyEnum'] },
                {
                    'tag': 'enumvalue',
                    'id': 'id_a',
                    'items': [
                        { 'tag': 'name', 'items': ['a'] },
                        {
                            'tag': 'briefdescription',
                            'items': ['an enumerator'],
                        },
                        {
                            'tag': 'detaileddescription',
                            'items': ['this is an enumerator'],
                        },
                    ]
                },
            ]
        }),
        None,
        None)
    e.resolve_references()
    assert len(e.objects) == 1
    e.objects[0].resolve_references()
    assert e.objects[0].id == 'id_a'
    assert e.objects[0].name == 'a'
    assert e.objects[0].enum == e
    assert e.objects[0].scope == e
    assert e.objects[0].type.text == e.name
    assert not e.objects[0].value.text
    assert e.objects[0].is_constexpr
    assert e.objects[0].is_const
    assert e.objects[0].is_static
    assert not e.objects[0].is_volatile
    assert not e.objects[0].is_inline
    assert collect_paragraphs(e.objects[0].brief) == 'an enumerator'
    assert (
        collect_paragraphs(e.objects[0].description)
        == 'this is an enumerator')
    assert len(e.members) == 1
    assert e.members[e.objects[0].name] == e.objects[0]

    e = docca.Enum(
        make_elem({
            'tag': 'compound',
            'id': 'enum1',
            'strong': 'yes',
            'items': [
                { 'tag': 'name', 'items': ['MyEnum'] },
                {
                    'tag': 'enumvalue',
                    'id': 'id_a',
                    'items': [
                        { 'tag': 'name', 'items': ['a'] },
                        {
                            'tag': 'briefdescription',
                            'items': ['an enumerator'],
                        },
                        {
                            'tag': 'detaileddescription',
                            'items': ['this is an enumerator'],
                        },
                    ]
                },
                {
                    'tag': 'enumvalue',
                    'id': 'id_b',
                    'items': [
                        { 'tag': 'name', 'items': ['b'] },
                        { 'tag': 'initializer', 'items': ['= 2'] },
                    ]
                },
            ]
        }),
        None,
        None)
    e.resolve_references()
    assert len(e.objects) == 2
    for o in e.objects:
        o.resolve_references()
    assert e.objects[1].id == 'id_b'
    assert e.objects[1].name == 'b'
    assert e.objects[1].value.text == '= 2'
    assert len(e.members) == 2
    assert e.members[e.objects[1].name] == e.objects[1]

    ns =  docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [ { 'tag': 'compoundname', 'items': ['ns'] } ]
        }))
    ns.resolve_references()

    e = docca.Enum(
        make_elem({
            'tag': 'compound',
            'id': 'enum1',
            'items': [
                { 'tag': 'name', 'items': ['MyEnum'] },
                {
                    'tag': 'enumvalue',
                    'id': 'id_a',
                    'items': [ { 'tag': 'name', 'items': ['a'] } ]
                },
                {
                    'tag': 'enumvalue',
                    'id': 'id_b',
                    'items': [ { 'tag': 'name', 'items': ['b'] } ]
                },
            ]
        }),
        None,
        ns)
    e.resolve_references()
    assert not e.is_scoped
    assert e.scope == ns
    assert e.objects[0].enum == e
    assert e.objects[0].scope == ns
    assert e.objects[1].enum == e
    assert e.objects[1].scope == ns

def test_type_alias():
    t = docca.TypeAlias(
        make_elem({
            'tag': 'memberdef',
            'id': 'alias1',
            'items': [
                { 'tag': 'name', 'items': ['MyAlias'] },
                { 'tag': 'type', 'items': ['std::string'] },
                {
                    'tag': 'briefdescription',
                    'items': ['a brief description'],
                },
                {
                    'tag': 'detaileddescription',
                    'items': ['a long description'],
                },
            ],
        }),
        None,
        None)
    t.resolve_references()
    assert t.declarator == 'using'
    assert t.id == 'alias1'
    assert t.name == 'MyAlias'
    assert t.aliased.text == 'std::string'
    assert collect_paragraphs(t.brief) == 'a brief description'
    assert collect_paragraphs(t.description) == 'a long description'

    index = {t.id: t}
    t2 = docca.TypeAlias(
        make_elem({
            'tag': 'memberdef',
            'id': 'alias2',
            'items': [
                { 'tag': 'name', 'items': ['MyAlias2'] },
                {
                    'tag': 'type',
                    'items': [ { 'tag': 'ref', 'refid': t.id } ],
                },
            ],
        }),
        None,
        None,
        index)
    t2.resolve_references()
    assert t2.id == 'alias2'
    assert t2.name == 'MyAlias2'
    assert t2.aliased.text == 'MyAlias'

def test_variable():
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'items': [
                                { 'tag': 'name', 'items': ['v1'] },
                                { 'tag': 'type', 'items': ['std::string'] },
                            ],
                        },
                        {
                            'tag': 'memberdef',
                            'id': 'idv2',
                            'kind': 'variable',
                            'items': [
                                { 'tag': 'name', 'items': ['v2'] },
                                { 'tag': 'type', 'items': ['std::string(*'] },
                                { 'tag': 'argsstring', 'items': [')(int)'] },
                            ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    assert len(ns.members) == 2
    for m in ns.members.values():
        m.resolve_references()
    assert isinstance(ns.members['v1'], docca.Variable)
    assert ns.members['v1'].id == 'idv1'
    assert ns.members['v1'].name == 'v1'
    assert ns.members['v1'].value.text == ''
    assert ns.members['v1'].type.text == 'std::string'
    assert not ns.members['v1'].is_static
    assert not ns.members['v1'].is_constexpr
    assert not ns.members['v1'].is_volatile
    assert not ns.members['v1'].is_const
    assert not ns.members['v1'].is_inline
    assert ns.members['v1'].access == 'public'

    assert isinstance(ns.members['v2'], docca.Variable)
    assert ns.members['v2'].id == 'idv2'
    assert ns.members['v2'].name == 'v2'
    assert ns.members['v2'].value.text == ''
    assert ns.members['v2'].type.text == 'std::string'
    assert ns.members['v2'].args.text == '(int)'
    assert not ns.members['v2'].is_static
    assert not ns.members['v2'].is_constexpr
    assert not ns.members['v2'].is_volatile
    assert not ns.members['v2'].is_const
    assert not ns.members['v2'].is_inline
    assert ns.members['v2'].access == 'public'

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'items': [
                                { 'tag': 'name', 'items': ['v1'] },
                                { 'tag': 'initializer', 'items': ['= 13'] },
                            ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].value.text == '= 13'

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'static': 'yes',
                            'items': [ { 'tag': 'name', 'items': ['v1'] } ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].is_static

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'inline': 'yes',
                            'items': [ { 'tag': 'name', 'items': ['v1'] } ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].is_inline

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'volatile': 'yes',
                            'items': [ { 'tag': 'name', 'items': ['v1'] } ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].is_volatile

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'const': 'yes',
                            'items': [ { 'tag': 'name', 'items': ['v1'] } ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].is_const

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'mutable': 'no',
                            'items': [ { 'tag': 'name', 'items': ['v1'] } ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].is_const

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'mutable': 'no',
                            'items': [ { 'tag': 'name', 'items': ['v1'] } ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].is_const

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'prot': 'private',
                            'items': [ { 'tag': 'name', 'items': ['v1'] } ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert ns.members['v1'].access == 'private'

    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns_id',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'idv1',
                            'kind': 'variable',
                            'items': [
                                { 'tag': 'name', 'items': ['v1'] },
                                {
                                    'tag': 'briefdescription',
                                    'items': ['!!!'],
                                },
                                {
                                    'tag': 'detaileddescription',
                                    'items': ['????'],
                                },
                            ],
                        },
                    ],
                },
            ]
        }))
    ns.resolve_references()
    for m in ns.members.values():
        m.resolve_references()
    assert collect_paragraphs(ns.members['v1'].brief) == '!!!'
    assert collect_paragraphs(ns.members['v1'].description) == '????'

def test_templatable():
    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [ { 'tag': 'compoundname', 'items': ['MyClass'] } ]
        }))
    c.resolve_references()
    assert not c.template_parameters
    assert not c.is_specialization

    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyClass'] },
                {
                    'tag': 'templateparamlist',
                    'items': [
                        {
                            'tag': 'param',
                            'items': [ { 'tag': 'type', 'items': ['class'] } ]
                        }
                    ]
                },
            ]
        }))
    c.resolve_references()
    assert len(c.template_parameters) == 1
    assert c.template_parameters[0].type.text == 'class'
    assert not c.template_parameters[0].description

    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyClass'] },
                {
                    'tag': 'templateparamlist',
                    'items': [
                        {
                            'tag': 'param',
                            'items': [
                                {
                                    'tag': 'defval',
                                    'items': ['= std::String']
                                }
                            ]
                        }
                    ]
                },
            ]
        }))
    c.resolve_references()
    assert c.template_parameters[0].default_value.text == '= std::String'

    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyClass'] },
                {
                    'tag': 'templateparamlist',
                    'items': [
                        {
                            'tag': 'param',
                            'items': [ { 'tag': 'declname', 'items': ['T'] } ]
                        }
                    ]
                },
            ]
        }))
    c.resolve_references()
    assert c.template_parameters[0].name == 'T'

    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyClass'] },
                {
                    'tag': 'templateparamlist',
                    'items': [
                        {
                            'tag': 'param',
                            'items': [
                                { 'tag': 'briefdescription', 'items': ['!'] },
                            ]
                        }
                    ]
                },
            ]
        }))
    c.resolve_references()
    assert (
        collect_paragraphs(c.template_parameters[0].description) == '!')

    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyClass'] },
                {
                    'tag': 'templateparamlist',
                    'items': [
                        {
                            'tag': 'param',
                            'items': [
                                { 'tag': 'type', 'items': ['std::size_t(&)'] },
                                { 'tag': 'array', 'items': ['[3]'] },
                            ]
                        }
                    ]
                },
            ]
        }))
    c.resolve_references()
    assert c.template_parameters[0].type.text == 'std::size_t'
    assert c.template_parameters[0].array.text == '[3]'

    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'someid',
            'items': [
                { 'tag': 'compoundname', 'items': ['MyClass<A, B, T>'] },
            ]
        }))
    c.resolve_references()
    assert c.is_specialization

def test_function():
    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'cl',
            'items': [ { 'tag': 'compoundname', 'items': ['Class'] } ]
        }))

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
            ]
        }),
        make_elem({}),
        c)
    c.members['f1'] = docca.OverloadSet([func])
    c.resolve_references()
    func.resolve_references()

    assert func.id == 'f1'
    assert func.name == 'func1'
    assert func.scope == c
    assert func.fully_qualified_name == 'Class::func1'
    assert not func.is_free
    assert not func.is_static
    assert not func.is_constexpr
    assert not func.is_volatile
    assert not func.is_const
    assert not func.is_friend
    assert not func.refqual
    assert not func.is_constructor
    assert not func.is_destructor
    assert not func.is_noexcept
    assert not func.is_explicit
    assert not func.is_specialization
    assert not func.parameters
    assert not func.template_parameters
    assert collect_paragraphs(func.brief) == ''
    assert collect_paragraphs(func.description) == ''
    assert func.access == 'public'
    assert func.virtual_kind == 'non-virtual'
    assert func.return_type.text == 'void'
    assert func.kind == 'nonstatic'

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'explicit': 'yes',
            'noexcept': 'yes',
            'static': 'yes',
            'constexpr': 'yes',
            'volatile': 'yes',
            'const': 'yes',
            'refqual': '&&',
            'virt': 'virtual',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['constexpr void'] },
                { 'tag': 'briefdescription', 'items': ['a function'] },
                { 'tag': 'detaileddescription', 'items': ['this function'] },
            ]
        }),
        make_elem({}),
        c)
    func.resolve_references()
    assert func.is_static
    assert func.is_constexpr
    assert func.is_volatile
    assert func.is_const
    assert func.refqual == '&&'
    assert func.is_noexcept
    assert func.is_explicit
    assert func.virtual_kind == 'virtual'
    assert func.kind == 'static'
    assert collect_paragraphs(func.brief) == 'a function'
    assert collect_paragraphs(func.description) == 'this function'

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
            ]
        }),
        make_elem({ 'kind': 'friend' }),
        c)
    func.resolve_references()
    assert func.is_friend
    assert func.kind == 'friend'

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
            ]
        }),
        make_elem({ 'kind': 'related' }),
        c)
    func.resolve_references()
    assert func.is_free
    assert func.kind == 'free'

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['() =delete'] },
                { 'tag': 'type', 'items': ['void'] },
            ]
        }),
        make_elem({}),
        c)
    func.resolve_references()
    assert func.is_deleted

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['() =default'] },
                { 'tag': 'type', 'items': ['void'] },
            ]
        }),
        make_elem({}),
        c)
    func.resolve_references()
    assert func.is_defaulted

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['Class'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
            ]
        }),
        make_elem({}),
        c)
    func.resolve_references()
    assert func.is_constructor

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['~Class'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
            ]
        }),
        make_elem({}),
        c)
    func.resolve_references()
    assert func.is_destructor
    assert func.is_noexcept

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                {
                    'tag': 'type',
                    'items': [ { 'tag': 'ref', 'refid': 'cl'} ]
                },
            ]
        }),
        make_elem({}),
        c,
        dict(cl=c))
    func.resolve_references()
    assert func.return_type.text == 'Class'

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
                {
                    'tag': 'param',
                    'items': [
                        { 'tag': 'type', 'items': ['int'] },
                        { 'tag': 'defval', 'items': ['= 12'] },
                        { 'tag': 'declname', 'items': ['A'] },
                        { 'tag': 'briefdescription', 'items': ['a param'] },
                    ],
                },
                {
                    'tag': 'param',
                    'items': [
                        {
                            'tag': 'type',
                            'items': [ { 'tag': 'ref', 'refid': 'cl'} ]
                        },
                    ],
                },
            ]
        }),
        make_elem({}),
        c,
        dict(cl=c))
    func.resolve_references()
    assert len(func.parameters) == 2

    assert func.parameters[0].type.text == 'int'
    assert func.parameters[0].default_value.text == '= 12'
    assert func.parameters[0].name == 'A'
    assert not func.parameters[0].array
    assert collect_paragraphs(func.parameters[0].description) == 'a param'

    assert func.parameters[1].type.text == 'Class'
    assert func.parameters[1].default_value.text == ''
    assert not func.parameters[1].name
    assert not func.parameters[1].array
    assert collect_paragraphs(func.parameters[1].description) == ''

    func = docca.Function(
        make_elem({
            'tag': 'memberdef',
            'id': 'f1',
            'items': [
                { 'tag': 'name', 'items': ['func1'] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
                {
                    'tag': 'param',
                    'items': [
                        { 'tag': 'type', 'items': ['int(&)'] },
                        { 'tag': 'array', 'items': ['[1]'] },
                    ],
                },
                {
                    'tag': 'param',
                    'items': [
                        { 'tag': 'type', 'items': ['int(*'] },
                        { 'tag': 'argsstring', 'items': [')(long)'] },
                    ],
                },
            ]
        }),
        make_elem({}),
        c)
    func.resolve_references()
    assert func.parameters[0].type.text == 'int'
    assert func.parameters[0].array.text == '[1]'
    assert func.parameters[1].type.text == 'int'
    assert func.parameters[1].args.text == '(long)'

    index = {}
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['Namespace'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'f1',
                            'kind': 'function',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': ['()'] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                },
            ],
        }),
        index);
    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'cl',
            'items': [
                { 'tag': 'compoundname', 'items': ['Class'] },
                {
                    'tag': 'sectiondef',
                    'kind': 'related',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'f1',
                            'kind': 'function',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': ['()'] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                },
            ],
        }),
        index);
    for e in index.values():
        e.resolve_references()
    assert len(c.members) == 1
    assert len(ns.members) == 0
    assert list(c.members.values())[0].scope == c

    index = {}
    c = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'cl',
            'items': [
                { 'tag': 'compoundname', 'items': ['Class'] },
                {
                    'tag': 'sectiondef',
                    'kind': 'related',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'f1',
                            'kind': 'function',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': ['()'] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                },
            ],
        }),
        index);
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['Namespace'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'id': 'f1',
                            'kind': 'function',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': ['()'] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                },
            ],
        }),
        index);
    for e in index.values():
        e.resolve_references()
    assert len(c.members) == 1
    assert len(ns.members) == 0
    assert list(c.members.values())[0].scope == c

def test_overload_set():
    index = dict()
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f1',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                }
            ]
        }))
    for entity in index.values():
        entity.resolve_references()
    assert len(ns.members) == 1
    oset = list(ns.members.values())[0]
    assert isinstance(oset, docca.OverloadSet)
    assert len(oset) == 1
    assert isinstance(oset[0], docca.Function)
    assert oset[0].overload_set == oset
    assert oset[0].is_sole_overload

    index = dict()
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f1',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f2',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                }
            ]
        }),
        index)
    for entity in index.values():
        entity.resolve_references()
    assert len(ns.members) == 1
    oset = list(ns.members.values())[0]
    assert len(oset) == 2
    assert oset[0].overload_set == oset
    assert not oset[0].is_sole_overload
    assert oset[0].overload_index == 0
    assert oset[1].overload_set == oset
    assert not oset[1].is_sole_overload
    assert oset[1].overload_index == 1

    index = dict()
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f1',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f2',
                            'items': [
                                { 'tag': 'name', 'items': ['func2'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                }
            ]
        }),
        index)
    for entity in index.values():
        entity.resolve_references()
    assert len(ns.members) == 2

    oset = list(ns.members.values())[0]
    assert len(oset) == 1
    assert oset[0].overload_set == oset
    assert oset[0].is_sole_overload

    oset = list(ns.members.values())[1]
    assert oset[0].overload_set == oset
    assert oset[0].is_sole_overload

    index = dict()
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f1',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f2',
                            'static': 'yes',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                            ],
                        },
                    ],
                }
            ]
        }),
        index)
    for entity in index.values():
        entity.resolve_references()
    oset = list(ns.members.values())[0]
    assert len(oset) == 1
    assert oset[0].overload_set == oset
    assert oset[0].is_sole_overload

    oset = list(ns.members.values())[1]
    assert oset[0].overload_set == oset
    assert oset[0].is_sole_overload

    index = dict()
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                {
                    'tag': 'sectiondef',
                    'items': [
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f1',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['void'] },
                                {
                                    'tag': 'briefdescription',
                                    'items': ['brief1'],
                                },
                            ],
                        },
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f2',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['int'] },
                                {
                                    'tag': 'briefdescription',
                                    'items': ['brief2'],
                                },
                            ],
                        },
                        {
                            'tag': 'memberdef',
                            'kind': 'function',
                            'id': 'f3',
                            'items': [
                                { 'tag': 'name', 'items': ['func1'] },
                                { 'tag': 'argsstring', 'items': [''] },
                                { 'tag': 'type', 'items': ['int'] },
                                {
                                    'tag': 'briefdescription',
                                    'items': ['brief1'],
                                },
                            ],
                        },
                    ],
                }
            ]
        }),
        index)
    for entity in index.values():
        entity.resolve_references()
    oset = list(ns.members.values())[0]
    assert len(oset) == 3
    assert collect_paragraphs(oset[0].brief) == 'brief1'
    assert oset[0].return_type.text == 'void'

    assert collect_paragraphs(oset[1].brief) == 'brief1'
    assert oset[1].return_type.text == 'int'

    assert collect_paragraphs(oset[2].brief) == 'brief2'

    funcs = sorted(list(oset))
    assert collect_paragraphs(funcs[0].brief) == 'brief1'
    assert funcs[0].return_type.text == 'void'

    assert collect_paragraphs(funcs[1].brief) == 'brief1'
    assert funcs[1].return_type.text == 'int'

    assert collect_paragraphs(funcs[2].brief) == 'brief2'

def test_parse_args():
    args = docca.parse_args([''])
    assert args.input is None
    assert args.output is None
    assert args.template is None
    assert args.directory is None
    assert len(args.config) == 0
    assert len(args.include) == 0
    assert len(args.extension) == 0

    with pytest.raises(SystemExit) as e:
        docca.parse_args(['', '-G', 'some string'])

    for flag in ('-i', '--input'):
        args = docca.parse_args(['', flag, 'some/input/file'])
        assert args.input == 'some/input/file'
        assert args.output is None
        assert args.template is None
        assert args.directory is None
        assert len(args.config) == 0
        assert len(args.include) == 0
        assert len(args.extension) == 0

        with pytest.raises(SystemExit) as e:
            docca.parse_args(
                ['', flag, 'some/input/file', flag, 'other/input/file'])

    for flag in ('-o', '--output'):
        args = docca.parse_args(['', flag, 'some/output/file'])
        assert args.input is None
        assert args.output == 'some/output/file'
        assert args.template is None
        assert args.directory is None
        assert len(args.config) == 0
        assert len(args.include) == 0
        assert len(args.extension) == 0

        with pytest.raises(SystemExit) as e:
            docca.parse_args(
                ['', flag, 'some/output/file', flag, 'other/output/file'])

    for flag in ('-T', '--template'):
        args = docca.parse_args(['', flag, 'some/template'])
        assert args.input is None
        assert args.output is None
        assert args.template  == 'some/template'
        assert args.directory is None
        assert len(args.config) == 0
        assert len(args.include) == 0
        assert len(args.extension) == 0

        with pytest.raises(SystemExit) as e:
            docca.parse_args(
                ['', flag, 'some/template', flag, 'other/template'])

    for flag in ('-D', '--directory'):
        args = docca.parse_args(['', flag, 'some/directory'])
        assert args.input is None
        assert args.output is None
        assert args.template is None
        assert args.directory == 'some/directory'
        assert len(args.config) == 0
        assert len(args.include) == 0
        assert len(args.extension) == 0

        with pytest.raises(SystemExit) as e:
            docca.parse_args(
                ['', flag, 'some/directory', flag, 'other/directory'])

    for flag in ('-c', '--config'):
        configs = []
        flags = []
        for n in range(5):
            configs.append( str(n) )
            flags.extend([ flag, configs[-1] ])
            args = docca.parse_args([''] + flags)
            assert args.input is None
            assert args.output is None
            assert args.template is None
            assert args.directory is None
            assert args.config == configs
            assert len(args.include) == 0
            assert len(args.extension) == 0

    for flag in ('-I', '--include'):
        includes = []
        flags = []
        for n in range(5):
            includes.append( str(n) )
            flags.extend([ flag, includes[-1] ])
            args = docca.parse_args([''] + flags)
            assert args.input is None
            assert args.output is None
            assert args.template is None
            assert args.directory is None
            assert len(args.config)  == 0
            assert args.include == includes
            assert len(args.extension) == 0

    for flag in ('-E', '--extension'):
        exts = []
        flags = []
        for n in range(5):
            exts.append( 'm' + str(n) )
            flags.extend([ flag, exts[-1] ])
            args = docca.parse_args([''] + flags)
            assert args.input is None
            assert args.output is None
            assert args.template is None
            assert args.directory is None
            assert len(args.config) == 0
            assert len(args.include) == 0
            assert args.extension == exts

    args = docca.parse_args([
        '', '-iinput1', '-ooutput2', '-Ttemplate3', '-cconf4', '-cconf5',
        '-Iinclude6', '-Iinclude7', '-Iinclude8', '-Ddir9', '-Eext10',
        '-Eext11',
    ])
    assert args.input == 'input1'
    assert args.output == 'output2'
    assert args.template == 'template3'
    assert args.directory == 'dir9'
    assert args.config == ['conf4', 'conf5']
    assert args.include == ['include6', 'include7', 'include8']
    assert args.extension == ['ext10', 'ext11']

def test_open_input(tmpdir):
    stdin = io.StringIO()
    cwd = '/no/such/path'

    args = argparse.Namespace()
    args.directory = None
    args.input = None

    file, ctx, dir = docca.open_input(stdin, args, cwd)
    assert file == stdin
    assert isinstance(ctx, docca.Nullcontext)
    assert dir == cwd

    args.directory = 'no/such/path'
    file, ctx, dir = docca.open_input(stdin, args, cwd)
    assert file == stdin
    assert isinstance(ctx, docca.Nullcontext)
    assert dir == 'no/such/path'

    text = '{ "key": "значение" }'
    args.directory = None
    args.input = os.path.join(tmpdir, 'input')
    with open(args.input, 'w', encoding='utf-8') as file:
        file.write(text)
    file, ctx, dir = docca.open_input(stdin, args, cwd)
    assert file.read() == text
    assert ctx == file
    assert dir == str(tmpdir)
    file.close()

    args.directory = 'no/such/path'
    file, ctx, dir = docca.open_input(stdin, args, cwd)
    assert file.read() == text
    assert ctx == file
    assert dir == 'no/such/path'
    file.close()

def test_collect_compond_refs():
    assert (
        list( docca.collect_compound_refs(io.StringIO(_compound_index)) )
        == [str(n) for n in range(5)])

def test_collect_data(tmpdir):
    kinds = [
        ('class', docca.Class,         'compoundname'),
        ('namespace', docca.Namespace, 'compoundname'),
        ('struct', docca.Struct,       'compoundname'),
        ('union', docca.Union,         'compoundname'),
        ('group', docca.Group,         'title'),
    ]
    for n, kind in enumerate(kinds):
        with open(os.path.join(tmpdir, str(n) + '.xml'), 'w') as f:
            f.write(_compound.format(n, kind[0], kind[2]))

    refs = list( docca.collect_compound_refs(io.StringIO(_compound_index)) )
    data = docca.collect_data(tmpdir, refs)
    for n, kind in enumerate(kinds):
        ref = str(n)
        assert isinstance(data[ref], kind[1])
        assert data[ref].id == ref
        assert data[ref].name == 'Entity' + ref

def test_open_output(tmpdir):
    stdout = io.StringIO()

    args = argparse.Namespace()
    args.output = None

    file, ctx = docca.open_output(stdout, args)
    assert file == stdout
    assert isinstance(ctx, docca.Nullcontext)

    args.output = os.path.join(tmpdir, 'output')
    file, ctx = docca.open_output(stdout, args)
    assert ctx == file
    with ctx:
        file.write('пример')

    with open(args.output, 'r', encoding='utf-8') as file:
        assert file.read() == 'пример'

def test_load_configs(tmpdir):
    args = argparse.Namespace()
    args.config = []

    conf = docca.load_configs(args)
    assert conf == {
        'include_private': False,
        'legacy_behavior': True,
    }

    c = os.path.join(tmpdir, 'c1')
    with open(c, 'w', encoding='utf-8') as f:
        f.write('{ "A": 1 }')
    args.config.append(c)
    conf = docca.load_configs(args)
    assert conf == {
        'include_private': False,
        'legacy_behavior': True,
        'A': 1,
    }

    c = os.path.join(tmpdir, 'c2')
    with open(c, 'w', encoding='utf-8') as f:
        f.write('{ "B": 2 }')
    args.config.append(c)
    conf = docca.load_configs(args)
    assert conf == {
        'include_private': False,
        'legacy_behavior': True,
        'A': 1,
        'B': 2,
    }

    c = os.path.join(tmpdir, 'c3')
    with open(c, 'w', encoding='utf-8') as f:
        f.write('{ "A": 3, "C": true }')
    args.config.append(c)
    conf = docca.load_configs(args)
    assert conf == {
        'include_private': False,
        'legacy_behavior': True,
        'A': 3,
        'B': 2,
        'C': True,
    }

def test_docca_include_dir():
    script = os.path.join('no/such/path/foo')
    assert docca.docca_include_dir(script) == 'no/such/path/include'

def test_template_file_name():
    args = argparse.Namespace()
    args.template = None

    assert (
        docca.template_file_name('no/such/path', args)
        == 'no/such/path/docca/quickbook.jinja2')

    args.template = 'no/such/path/too/tmpl'
    assert docca.template_file_name('no/such/path', args) == args.template

def test_collect_include_dirs():
    args = argparse.Namespace()
    args.include = []

    assert (
        docca.collect_include_dirs('path/to/tmpl', 'docca/include/dir', args)
        == ['path/to', 'docca/include/dir'])

    args.include = ['1/2', '3/4']
    assert (
        docca.collect_include_dirs('path/to/tmpl', 'docca/include/dir', args)
        == ['path/to', 'docca/include/dir', '1/2', '3/4'])

def test_construct_environment():
    conf = dict()
    loader = jinja2.DictLoader({})
    env = docca.construct_environment(loader, conf)

    assert not env.autoescape
    assert env.loader == loader
    assert env.undefined == jinja2.StrictUndefined

    assert [
        ext for ext in env.extensions.values()
        if isinstance(ext, jinja2.ext.do)
    ]
    assert [
        ext for ext in env.extensions.values()
        if isinstance(ext, jinja2.ext.loopcontrols)
    ]

    assert env.globals['Access'] == docca.Access
    assert env.globals['FunctionKind'] == docca.FunctionKind
    assert env.globals['VirtualKind'] == docca.VirtualKind
    assert env.globals['Section'] == docca.Section
    assert env.globals['ParameterList'] == docca.ParameterList
    assert env.globals['Config'] == conf
    assert env.globals['re'] == re

def test_load_extensions():
    exts = docca.load_extensions([])
    assert len(exts) == 0

    here = os.path.dirname(__file__)
    exts = docca.load_extensions([ os.path.join(here, 'exts/1.py') ])
    assert len(exts) == 1

    assert exts[0].__file__ == os.path.join(here, 'exts/1.py')
    assert exts[0].__name__ == 'docca._ext0'

    exts = docca.load_extensions([
        os.path.join(here, 'exts/2.py'), os.path.join(here, 'exts/3.py') ])
    assert len(exts) == 2
    assert exts[0].__file__ == os.path.join(here, 'exts/2.py')
    assert exts[0].__name__ == 'docca._ext0'
    assert exts[1].__file__ == os.path.join(here, 'exts/3.py')
    assert exts[1].__name__ == 'docca._ext1'

def test_install_extensions():
    env = docca.construct_environment(dict(), jinja2.DictLoader({}))

    default_globals = env.globals.copy()
    default_tests = env.tests.copy()
    default_filters = env.filters.copy()
    default_extensions = env.extensions.copy()

    env = docca.install_extensions(env, [])
    assert env.globals == default_globals
    assert env.filters == default_filters
    assert env.tests == default_tests
    assert env.extensions == default_extensions

    class mod1:
        @staticmethod
        def install_docca_extension(env):
            env.globals['foo'] = 'bar'
    env = docca.construct_environment(dict(), jinja2.DictLoader({}))
    env = docca.install_extensions(env, [mod1])
    assert env.globals['foo'] == 'bar'

    d = lambda x: False

    class mod2:
        @staticmethod
        def install_docca_extension(env):
            env.globals['a'] = 'b'
            env.tests['c'] = d
    env = docca.construct_environment(dict(), jinja2.DictLoader({}))
    env = docca.install_extensions(env, [mod1, mod2])
    assert env.globals['foo'] == 'bar'
    assert env.globals['a'] == 'b'
    assert env.tests['c'] == d

def test_render():
    file = io.StringIO()
    loader = jinja2.DictLoader({'tmpl': _simple_template})
    env = jinja2.Environment(loader=loader)
    docca.render(env, 'path/to/tmpl', file, [1,2,3,4])
    assert file.getvalue() == '1, 2, 3, 4'

_compound_index = '''\
<?xml version='1.0'?>
<doxygenindex>
  <compound refid="0" kind="class"/>
  <compound refid="B" kind="dir"/>
  <compound refid="1" kind="class"/>
  <compound refid="2" kind="class"/>
  <compound refid="3" kind="class"/>
  <compound refid="A" kind="file"/>
  <compound refid="4" kind="namespace"/>
</doxygenindex>
'''
_compound = '''\
<?xml version='1.0'?>
<doxygen>
  <compounddef id="{0}" kind="{1}" prot="public">
    <{2}>Entity{0}</{2}>
  </compounddef>
</doxygen>
'''
_simple_template = '''\
{%- set sep = joiner(", ") -%}
{% for item in entities %}{{sep()}}{{item}}{% endfor %}
'''
