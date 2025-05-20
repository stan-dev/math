#
# Copyright (c) 2024 Dmitry Arkhipov (grisumbras@yandex.ru)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

import docca
import io
import jinja2
import os
import pytest
import textwrap

from docca_test_helpers import make_elem


@pytest.fixture
def cfg():
    return dict(legacy_behavior=False)

@pytest.fixture
def entities():
    index = dict()

    f1 = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'f1',
        'items': [
            { 'tag': 'name', 'items': ['func'] },
            { 'tag': 'argsstring', 'items': ['()'] },
            { 'tag': 'type', 'items': ['void'] },
            {
                'tag': 'param',
                'items': [ { 'tag': 'type', 'items': ['int'] } ]
            },
            {
                'tag': 'param',
                'items': [
                    { 'tag': 'type', 'items': ['int(&)'] },
                    { 'tag': 'declname', 'items': ['n'] },
                    { 'tag': 'array', 'items': ['[1]'] },
                ],
            },
            {
                'tag': 'param',
                'items': [
                    {
                        'tag': 'type',
                        'items': [
                            { 'tag': 'ref', 'refid': 'cl1', 'items': ['cl2'] },
                        ],
                    },
                    {
                        'tag': 'defval',
                        'items': [
                            { 'tag': 'ref', 'refid': 'v1', 'items': ['v1'] },
                        ],
                    },
                ],
            },
            {
                'tag': 'param',
                'items': [
                    { 'tag': 'type', 'items': ['unsigned'] },
                    { 'tag': 'declname', 'items': ['t'] },
                    { 'tag': 'defval', 'items': ['2'] },
                ],
            },
        ]
    }
    f2 = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'f2',
        'static': 'yes',
        'constexpr': 'yes',
        'noexcept': 'yes',
        'items': [
            { 'tag': 'name', 'items': ['func'] },
            { 'tag': 'type', 'items': ['constexpr void'] },
            {
                'tag': 'argsstring',
                 'items': ['() noexcept(noexcept((void)12))'],
            },
        ]
    }
    g1 = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'g1',
        'refqual': 'rvalue',
        'noexcept': 'yes',
        'noexceptexpression': 'foobar(")", "(")',
        'items': [
            { 'tag': 'name', 'items': ['g'] },
            { 'tag': 'argsstring', 'items': ['()'] },
            { 'tag': 'type', 'items': ['void'] },
            { 'tag': 'briefdescription', 'items': ['One overload of g'] },
        ]
    }
    g2 = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'g2',
        'const': 'yes',
        'refqual': 'lvalue',
        'items': [
            { 'tag': 'name', 'items': ['g'] },
            { 'tag': 'argsstring', 'items': ['()'] },
            { 'tag': 'type', 'items': ['void'] },
            { 'tag': 'briefdescription', 'items': ['Another overload of g'] },
        ]
    }
    g3 = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'g3',
        'refqual': 'lvalue',
        'items': [
            { 'tag': 'name', 'items': ['g'] },
            { 'tag': 'argsstring', 'items': ['() =delete'] },
            { 'tag': 'type', 'items': ['void'] },
            { 'tag': 'briefdescription', 'items': ['Third overload of g'] },
        ]
    }

    cons = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'cl1_c',
        'explicit': 'yes',
        'constexpr': 'yes',
        'items': [
            { 'tag': 'name', 'items': ['klass'] },
            { 'tag': 'argsstring', 'items': ['()'] },
        ]
    }
    destr = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'cl1_d',
        'noexcept': 'no',
        'items': [
            { 'tag': 'name', 'items': ['~klass'] },
            { 'tag': 'argsstring', 'items': ['()'] },
        ]
    }

    def operator(s):
        return {
            'tag': 'memberdef',
            'kind': 'function',
            'id': 'o' + s,
            'noexcept': 'yes',
            'items': [
                { 'tag': 'name', 'items': ['operator' + s] },
                { 'tag': 'argsstring', 'items': ['()'] },
                { 'tag': 'type', 'items': ['void'] },
                {
                    'tag': 'briefdescription',
                    'items': [
                        'This is ',
                        { 'tag': 'emphasis', 'items': ['one'] },
                        ' of the operators',
                    ],
                },
            ]
        }

    cl4 = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'cl4',
            'items': [
                { 'tag': 'compoundname', 'items': ['nested2'] },
                {
                    'tag': 'detaileddescription',
                    'items': [
                        'See ',
                        { 'tag': 'ref', 'refid': 'cl4', 'items': ['nested2'] },
                        '.',
                    ],
                },
            ]
        }),
        index)
    cl2 = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'cl2',
            'items': [
                { 'tag': 'compoundname', 'items': ['nested'] },
                { 'tag': 'innerclass', 'refid': cl4.id },
                {
                    'tag': 'templateparamlist',
                    'items': [
                        {
                            'tag': 'param',
                            'items': [ { 'tag': 'type', 'items': ['class'] } ]
                        },
                        {
                            'tag': 'param',
                            'items': [
                                { 'tag': 'type', 'items': ['int(&)'] },
                                { 'tag': 'declname', 'items': ['N'] },
                                { 'tag': 'array', 'items': ['[1]'] },
                            ],
                        },
                        {
                            'tag': 'param',
                            'items': [
                                { 'tag': 'type', 'items': ['class'] },
                                { 'tag': 'defval', 'items': ['std::String'] },
                            ],
                        },
                        {
                            'tag': 'param',
                            'items': [
                                { 'tag': 'type', 'items': ['class'] },
                                { 'tag': 'declname', 'items': ['T'] },
                                { 'tag': 'defval', 'items': ['std::String'] },
                            ],
                        },
                    ]
                },
            ]
        }),
        index)
    cl1 = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'cl1',
            'items': [
                { 'tag': 'compoundname', 'items': ['klass'] },
                {
                    'tag': 'detaileddescription',
                    'items': [
                        'This class contains ',
                        { 'tag': 'ref', 'refid': 'cl2', 'items': ['nested'] },
                        '.',
                    ],
                },
                { 'tag': 'innerclass', 'refid': cl2.id },
                { 'tag': 'sectiondef', 'kind': 'friend', 'items': [f1] },
                {
                    'tag': 'sectiondef',
                    'kind': 'public',
                    'items': [f2, g1, g2, g3],
                },
                {
                    'tag': 'sectiondef', 'kind': 'public',
                    'items': [cons,destr],
                },
                {
                    'tag': 'sectiondef',
                    'kind': 'ops',
                    'items': [
                        operator(s)
                        for s in (
                            '=', '+=', '-=', '*=', '/=', '%=', '&=', '|=',
                            '^=', '<<=', '>>=', '++', '--', '+', '-', '*', '/',
                            '%', '~', '&', '|', '^', '<<', '>>', '!', '&&',
                            '||', '==', '!=', '<=', '>=', '>', '<', '<=>',
                            '[]', '->', '()', ','
                        )
                    ],
                },
            ]
        }),
        index)
    cl3 = docca.Class(
        make_elem({
            'tag': 'compound',
            'id': 'cl3',
            'items': [
                { 'tag': 'compoundname', 'items': ['tmpl<int, std::string>'] },
                {
                    'tag': 'basecompoundref',
                    'virt': 'virtual',
                    'prot': 'public',
                    'refid': 'cl2',
                },
            ]
        }),
        index)
    ns2 = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns2',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns2'] },
                { 'tag': 'innerclass', 'refid': cl1.id },
                { 'tag': 'innerclass', 'refid': cl3.id },
            ]
        }),
        index)

    e1 = {
        'tag': 'memberdef',
        'kind': 'enum',
        'id': 'e1',
        'strong': 'yes',
        'items': [
            { 'tag': 'name', 'items': ['e1'] },
            { 'tag': 'type', 'items': ['unsigned'] },
            {
                'tag': 'enumvalue',
                'id': 'e1_a',
                'items': [ { 'tag': 'name', 'items': ['a'] } ],
            },
        ]
    }
    t1 = {
        'tag': 'memberdef',
        'kind': 'typedef',
        'id': 't1',
        'items': [
            { 'tag': 'name', 'items': ['TA'] },
            { 'tag': 'type', 'items': ['std::string'] },
        ]
    }
    v1 = {
        'tag': 'memberdef',
        'kind': 'variable',
        'id': 'v1',
        'static': 'yes',
        'items': [
            { 'tag': 'name', 'items': ['Var'] },
            { 'tag': 'type', 'items': ['int'] },
            { 'tag': 'initializer', 'items': ['= 11'] }
        ]
    }
    v2 = {
        'tag': 'memberdef',
        'kind': 'variable',
        'id': 'v2',
        'static': 'yes',
        'constexpr': 'yes',
        'items': [
            { 'tag': 'name', 'items': ['Var2'] },
            { 'tag': 'type', 'items': ['double'] },
        ]
    }
    v3 = {
        'tag': 'memberdef',
        'kind': 'variable',
        'id': 'v3',
        'items': [
            { 'tag': 'name', 'items': ['Var3'] },
            { 'tag': 'type', 'items': ['int(*'] },
            { 'tag': 'argsstring', 'items': [')(int, long)'] },
        ]
    }

    ns1 = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns1',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns1'] },
                { 'tag': 'innernamespace', 'refid': ns2.id },
                { 'tag': 'sectiondef', 'items': [e1, t1, v1, v2, v3] },
                {
                    'tag': 'briefdescription',
                    'items': [
                        'This namespace contains ',
                        { 'tag': 'ref', 'refid': 'cl1', 'items': ['nested'] },
                        '.',
                    ],
                },
                {
                    'tag': 'location',
                    'file': '/path/to/ns1.hpp',
                    'line': '223',
                    'column': '322'
                },
            ]
        }),
        index)

    ostream = {
        'tag': 'memberdef',
        'kind': 'typedef',
        'id': 'ostream',
        'items': [
            { 'tag': 'name', 'items': ['ostream'] },
            { 'tag': 'type', 'items': ['void'] },
            { 'tag': 'briefdescription', 'items': ['!!!'] },
            {
                'tag': 'detaileddescription',
                'items': [
                    {
                        'tag': 'simplesect',
                        'kind': 'see',
                        'items': [
                            { 'tag': 'para', 'items': ['  http://ostream.org  '] },
                        ],
                    },
                ],
            }
        ],
    }
    std = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'std',
            'items': [
                { 'tag': 'compoundname', 'items': ['std'] },
                { 'tag': 'sectiondef', 'items': [ostream] },
            ]
        }),
        index)
    for e in index.values():
        e.resolve_references()

    return index


class Renderer():
    def __init__(self, cfg):
        self._templates = dict()
        self.loader = jinja2.ChoiceLoader([
            jinja2.DictLoader(self._templates),
            jinja2.FileSystemLoader(
                [os.path.join(os.path.dirname(__file__), '..', 'include')]),
        ])
        self.env = docca.construct_environment(self.loader, cfg)
        self.template = ''

    def __call__(self, data):
        file = io.StringIO()
        docca.render(self.env, '__tmpl__', file, data)
        return file.getvalue()

    @property
    def template(self):
        return self._templates['__tmpl__']

    @template.setter
    def template(self, tmpl):
        self._templates['__tmpl__'] = tmpl

@pytest.fixture
def render(cfg):
    return Renderer(cfg)


def test_escape(render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.escape(entities) }}'''
    assert render('[[foobar]]') == '\\[\\[foobar\\]\\]'

def test_text_helper(cfg, render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.text_helper(entities) }}'''
    assert render('') == ''

    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.text_helper(entities) }}'''
    assert render('[[foobar]]') == '\\[\\[foobar\\]\\]'

    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.text_helper(entities, in_code=True) }}'''
    cfg['replace_strings'] = {'foobar': 'FooBar'}
    assert render('[[foobar]]') == '[[FooBar]]'

    cfg['replace_strings'] = {'foob\\b': 'FooBar'}
    assert render('[[foobar]]') == '[[foobar]]'

    cfg['replace_strings'] = {'\\bf(o+)bar\\b': 'F\\1BaR'}
    assert render('[[fobar]]') == '[[FoBaR]]'
    assert render('[[foobar]]') == '[[FooBaR]]'
    assert render('[[fooobar]]') == '[[FoooBaR]]'

    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.text_helper(entities, in_code=True) }}'''
    assert render('const int   &')  == 'const int&'
    assert render('const int   &&') == 'const int&&'
    assert render('const int   *')  == 'const int*'
    assert render('const int   &...')  == 'const int&...'
    assert render('const int   &&...') == 'const int&&...'
    assert render('const int   *...')  == 'const int*...'

def test_link(entities, cfg, render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.link(entities) }}'''
    assert render(entities['ns1']) == 'ns1'
    assert render(entities['cl1']) == 'ns1__ns2__klass'
    assert render(entities['cl2']) == 'ns1__ns2__klass.nested'
    assert render(entities['cl3']) == 'ns1__ns2__tmpl_lt_int_comma__std__string_gt_'
    assert render(entities['e1']) == 'ns1__e1'
    assert render(entities['e1_a']) == 'ns1__e1'
    assert render(entities['f1']) == 'ns1__ns2__klass.func_fr'
    assert render(entities['f2']) == 'ns1__ns2__klass.func'
    assert render(entities['g1']) == 'ns1__ns2__klass.g.overload1'
    assert render(entities['g2']) == 'ns1__ns2__klass.g.overload2'
    assert render(entities['o=']) == 'ns1__ns2__klass.operator_eq_'
    assert render(entities['o+=']) == 'ns1__ns2__klass.operator_plus__eq_'
    assert render(entities['o-=']) == 'ns1__ns2__klass.operator_minus__eq_'
    assert render(entities['o*=']) == 'ns1__ns2__klass.operator__star__eq_'
    assert render(entities['o/=']) == 'ns1__ns2__klass.operator_slash__eq_'
    assert render(entities['o%=']) == 'ns1__ns2__klass.operator_mod__eq_'
    assert render(entities['o&=']) == 'ns1__ns2__klass.operator_and__eq_'
    assert render(entities['o|=']) == 'ns1__ns2__klass.operator_or__eq_'
    assert render(entities['o^=']) == 'ns1__ns2__klass.operator_xor__eq_'
    assert render(entities['o<<=']) == 'ns1__ns2__klass.operator_lt__lt__eq_'
    assert render(entities['o>>=']) == 'ns1__ns2__klass.operator__gt__gt__eq_'
    assert render(entities['o++']) == 'ns1__ns2__klass.operator_plus__plus_'
    assert render(entities['o--']) == 'ns1__ns2__klass.operator_minus__minus_'
    assert render(entities['o+']) == 'ns1__ns2__klass.operator_plus_'
    assert render(entities['o*']) == 'ns1__ns2__klass.operator__star_'
    assert render(entities['o/']) == 'ns1__ns2__klass.operator_slash_'
    assert render(entities['o%']) == 'ns1__ns2__klass.operator_mod_'
    assert render(entities['o~']) == 'ns1__ns2__klass.operator_bnot_'
    assert render(entities['o&']) == 'ns1__ns2__klass.operator_and_'
    assert render(entities['o|']) == 'ns1__ns2__klass.operator_or_'
    assert render(entities['o^']) == 'ns1__ns2__klass.operator_xor_'
    assert render(entities['o<<']) == 'ns1__ns2__klass.operator_lt__lt_'
    assert render(entities['o>>']) == 'ns1__ns2__klass.operator__gt__gt_'
    assert render(entities['o!']) == 'ns1__ns2__klass.operator__not_'
    assert render(entities['o&&']) == 'ns1__ns2__klass.operator_and__and_'
    assert render(entities['o||']) == 'ns1__ns2__klass.operator_or__or_'
    assert render(entities['o==']) == 'ns1__ns2__klass.operator_eq__eq_'
    assert render(entities['o!=']) == 'ns1__ns2__klass.operator__not__eq_'
    assert render(entities['o<=']) == 'ns1__ns2__klass.operator_lt__eq_'
    assert render(entities['o>=']) == 'ns1__ns2__klass.operator__gt__eq_'
    assert render(entities['o<']) == 'ns1__ns2__klass.operator_lt_'
    assert render(entities['o<=>']) == 'ns1__ns2__klass.operator_spshp_'
    assert render(entities['o[]']) == 'ns1__ns2__klass.operator__lb__rb_'
    assert render(entities['o->']) == 'ns1__ns2__klass.operator__arrow_'
    assert render(entities['o()']) == 'ns1__ns2__klass.operator_lp__rp_'
    assert render(entities['o,']) == 'ns1__ns2__klass.operator_comma_'

    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.link(entities, prefer_overload=True) }}'''
    assert render(entities['g1']) == 'ns1__ns2__klass.g'
    assert render(entities['g2']) == 'ns1__ns2__klass.g'

    cfg['link_prefix'] = 'a.b.c.'
    assert render(entities['cl1']) == 'a.b.c.ns1__ns2__klass'

def test_link(entities, render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.anchor(entities) }}'''
    assert render(entities['ns1']) == 'ns1'
    assert render(entities['cl1']) == 'ns1__ns2__klass'
    assert render(entities['cl2']) == 'nested'
    assert render(entities['cl3']) == 'ns1__ns2__tmpl_lt_int_comma__std__string_gt_'
    assert render(entities['e1']) == 'ns1__e1'
    assert render(entities['f1']) == 'func_fr'
    assert render(entities['f2']) == 'func'
    assert render(entities['g1']) == 'overload1'
    assert render(entities['g2']) == 'overload2'
    assert render(entities['o=']) == 'operator_eq_'

def test_abridged_fqn(entities, cfg, render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.abridged_fqn(entities) }}'''
    assert render(entities['ns1']) == 'ns1'
    assert render(entities['cl1']) == 'ns1::ns2::klass'
    assert render(entities['g1']) == 'ns1::ns2::klass::g'
    assert render(entities['o[]']) == 'ns1::ns2::klass::operator\\[\\]'

    cfg['default_namespace'] = 'ns3'
    assert render(entities['ns1']) == 'ns1'
    assert render(entities['cl1']) == 'ns1::ns2::klass'

    cfg['default_namespace'] = 'ns1'
    assert render(entities['ns1']) == 'ns1'
    assert render(entities['cl1']) == 'ns2::klass'

def test_phrase(entities, cfg, render):
    cfg['external_marker'] = '!!!'
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.phrase(entities) }}'''
    parts = [
        'regular text ',
        docca.Monospaced(['monospaced text']),
        ' ',
        docca.Emphasised(['emph ', docca.Strong(['bold'])]),
        docca.Linebreak(),
        docca.UrlLink('http://a.b/c', ['link text']),
        docca.EntityRef(entities['g1'], ['ref text']),
        docca.Linebreak(),
        docca.EntityRef(entities['ostream'], ['output stream']),
        docca.Phrase([' ', 'more text']),
        docca.Linebreak(),
        docca.EntityRef(entities['o[]'], 'operator[]'),
    ]
    assert render(parts) == textwrap.dedent('''\
        regular text `monospaced text` ['emph [*bold]]

        [@http://a.b/c link text][link ns1__ns2__klass.g `ref text`]

        [@http://ostream.org `output stream`] more text

        [link ns1__ns2__klass.operator__lb__rb_ `operator[]`]''')

    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.phrase(entities, in_code=True) }}'''
    assert render(parts) == textwrap.dedent('''\
        regular text `monospaced text` ['emph [*bold]]

        [@http://a.b/c link text]``[link ns1__ns2__klass.g [^ns1::ns2::klass::g]]``

        ``[@http://ostream.org [^std::ostream]]`` more text

        ``[link ns1__ns2__klass.operator__lb__rb_ [^ns1::ns2::klass::operator\[\]]]``''')

def test_description(render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.description(entities, nesting='  ') }}'''
    assert render([]) == ''

    parts = [
        docca.Paragraph([ '1 ', docca.Monospaced(['2']) ]),
        docca.CodeBlock([ 'int i = 0;', '++i;' ]),
        docca.List(
            '1',
            [
                docca.ListItem([
                    docca.Paragraph(['o1']),
                    docca.List(
                        None,
                        [
                            docca.ListItem([ docca.Paragraph(['u1']) ]),
                            docca.ListItem([ docca.Paragraph(['u2']) ]),
                        ]),
                ]),
                docca.ListItem([ docca.Paragraph(['o2']) ]),
            ]),
        docca.Section(
            None,
            docca.Paragraph([ 'Title' ]),
            [
                docca.Paragraph(['3']),
                docca.Paragraph(['4']),
            ]),
        docca.Table(
            None,
            [
                [
                    docca.Cell([ docca.Paragraph(['tl']) ]),
                    docca.Cell([ docca.Paragraph(['tr']) ]),
                ],
                [
                    docca.Cell([ docca.Paragraph(['bl']) ]),
                    docca.Cell([ docca.Paragraph(['br']) ]),
                ],
            ],
            caption=docca.Paragraph([ 'Caption' ])),
        docca.ParameterList(
            'param',
            [
                docca.ParameterDescription(
                    [docca.Paragraph(['param1 description'])],
                    [
                        docca.ParameterItem(
                            docca.Phrase(['int']),
                            docca.Phrase(['a']),
                            None),
                        docca.ParameterItem(
                            docca.Phrase(['char const*']),
                            docca.Phrase(['b']),
                            None),
                    ],
                ),
                docca.ParameterDescription(
                    [docca.Paragraph(['param2 description'])],
                    [
                        docca.ParameterItem(
                            docca.Phrase(['float']),
                            docca.Phrase(['c']),
                            None),
                    ]
                ),
            ]),
    ]
    assert render(parts) == textwrap.dedent('''\
        1 `2`


          ```
          int i = 0;
          ++i;
          ```
          # o1

                 * u1

                 * u2

          # o2

        [heading Title]
        3

        4

        [table Caption
          [
            [ tl

            ]
            [ tr

            ]
          ]
          [
            [ bl

            ]
            [ br

            ]
          ]
        ]
        [heading Parameters]
        [table [[Name][Description]]
          [
            [`int a`, `char const* b`
            ]
            [
        param1 description


            ]
          ]
          [
            [`float c`
            ]
            [
        param2 description


            ]
          ]
        ]

        ''')

    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.description(entities, nesting='  ', title='Title') }}'''
    assert render([]) == ''

    parts = [ docca.Paragraph([ '1 ', docca.Monospaced(['2']) ]) ]
    assert render(parts) == textwrap.dedent('''\
        [heading Title]
        1 `2`

        ''')
    parts = [
        docca.ParameterList(
            'param',
            [
                docca.ParameterDescription(
                    [docca.Paragraph(['param1 description'])],
                    [
                        docca.ParameterItem(
                            docca.Phrase(['int']),
                            docca.Phrase(['a']),
                            None),
                        docca.ParameterItem(
                            docca.Phrase(['char const*']),
                            docca.Phrase(['b']),
                            None),
                    ],
                ),
                docca.ParameterDescription(
                    [docca.Paragraph(['param2 description'])],
                    [
                        docca.ParameterItem(
                            docca.Phrase(['float']),
                            docca.Phrase(['c']),
                            None),
                    ]
                ),
            ])
    ]
    assert render(parts) == textwrap.dedent('''\
        [heading Parameters]
        [table [[Name][Description]]
          [
            [`int a`, `char const* b`
            ]
            [
        param1 description


            ]
          ]
          [
            [`float c`
            ]
            [
        param2 description


            ]
          ]
        ]

        ''')

def test_subsection(render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.subsection(entities) }}'''
    assert render(docca.Section('see', None, [])) == '[heading See Also]\n'
    assert render(docca.Section('return', None, [])) == '[heading Return Value]\n'
    assert render(docca.Section('author', None, [])) == '[heading Author]\n'
    assert render(docca.Section('authors', None, [])) == '[heading Authors]\n'
    assert render(docca.Section('version', None, [])) == '[heading Version]\n'
    assert render(docca.Section('since', None, [])) == '[heading Since]\n'
    assert render(docca.Section('date', None, [])) == '[heading Date]\n'
    assert render(docca.Section('note', None, [])) == '[heading Remarks]\n'
    assert render(docca.Section('warning', None, [])) == '[heading Warning]\n'
    assert render(docca.Section('pre', None, [])) == '[heading Preconditions]\n'
    assert render(docca.Section('post', None, [])) == '[heading Postconditions]\n'
    assert render(docca.Section('copyright', None, [])) == '[heading Copyright]\n'
    assert render(docca.Section('invariant', None, [])) == '[heading Invariants]\n'
    assert render(docca.Section('remark', None, [])) == '[heading Remarks]\n'
    assert render(docca.Section('attention', None, [])) == '[heading Attention]\n'
    assert render(docca.Section('par', None, [])) == '[heading Paragraph]\n'
    assert render(docca.Section('rcs', None, [])) == '[heading RCS]\n'

def test_parameter_list(render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.parameter_list(entities) }}'''
    params = docca.ParameterList(
        'param',
        [ docca.ParameterDescription(docca.Paragraph([]), []) ])
    assert render(params) == textwrap.dedent('''\
        [heading Parameters]
        [table [[Name][Description]]
          [
            [
            ]
            [

            ]
          ]
        ]
        ''')
    params.kind = 'retval'
    assert render(params) == textwrap.dedent('''\
        [heading Return Values]
        [table [[Type][Description]]
          [
            [
            ]
            [

            ]
          ]
        ]
        ''')
    params.kind = 'exception'
    assert render(params) == textwrap.dedent('''\
        [heading Exceptions]
        [table [[Type][Thrown On]]
          [
            [
            ]
            [

            ]
          ]
        ]
        ''')
    params.kind = 'templateparam'
    assert render(params) == textwrap.dedent('''\
        [heading Template Parameters]
        [table [[Type][Description]]
          [
            [
            ]
            [

            ]
          ]
        ]
        ''')

def test_heading(render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.heading(entities) }}'''
    assert render('some header') == '[heading some header]'

def test_summary_table(render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {% call(row) qbk.summary_table(entities[2:], entities[0], entities[1]) %}
        [-{{row}}-]
        {% endcall %}'''
    assert (
        render(['A', ['B'], 'C'])
        == '[heading A]\n'
           '[table [[B]]\n'
           '  [\n'
           '    \n'
           '        [-C-]\n'
           '        \n'
           '  ]\n'
           ']\n')
    assert (
        render(['TTT Tttt', ['A', 'B'], 'C', 'D'])
        == '[heading TTT Tttt]\n'
           '[table [[A][B]]\n'
           '  [\n'
           '    \n'
           '        [-C-]\n'
           '        \n'
           '  ]\n'
           '  [\n'
           '    \n'
           '        [-D-]\n'
           '        \n'
           '  ]\n'
           ']\n')

def test_simple_summary_table(entities, render):
    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.simple_summary_table(entities, 'T') }}'''
    assert render([entities['ns1'], entities['o[]']]) == textwrap.dedent('''\
        [heading T]
        [table [[Name][Description]]
          [
            [[*[link ns1 ns1]]
            ]
            [This namespace contains [link ns1__ns2__klass `nested`].
            ]
          ]
          [
            [[*[link ns1__ns2__klass.operator__lb__rb_ operator\\[\\]]]
            ]
            [This is ['one] of the operators
            ]
          ]
        ]
        ''')

def test_function_summary_table(entities, render):
    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.function_summary_table(entities, 'T') }}'''
    assert render([
            entities['cl1_c'].overload_set,
            entities['cl1_d'].overload_set,
            entities['o[]'].overload_set,
            entities['g1'].overload_set,
        ]) == textwrap.dedent('''\
            [heading T]
            [table [[Name][Description]]
              [
                [[*[link ns1__ns2__klass.g g]]
                ]
                [One overload of g \'''<sbr/>\'''[role silver \\u2014]\'''<sbr/>\'''Another overload of g \'''<sbr/>\'''[role silver \\u2014]\'''<sbr/>\'''Third overload of g
                ]
              ]
              [
                [[*[link ns1__ns2__klass.klass klass]\\u00A0[role silver \\[constructor\\]]]
                ]
                [
                ]
              ]
              [
                [[*[link ns1__ns2__klass.operator__lb__rb_ operator\\[\\]]]
                ]
                [This is [\'one] of the operators
                ]
              ]
              [
                [[*[link ns1__ns2__klass._dtor_klass ~klass]\\u00A0[role silver \\[destructor\\]]]
                ]
                [
                ]
              ]
            ]
            ''')

def test_template_parameters(entities, render):
    render.template = '''
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.template_parameters(entities) }}'''
    assert render(entities['cl1']) == ''
    assert render(entities['cl3']) == 'template<>\n'
    assert render(entities['cl2']) == textwrap.dedent('''\
        template<
            class,
            int (&N) [1],
            class = std::String,
            class T = std::String>
        ''')

def test_source_header(render):
    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.source_header(entities) }}'''
    assert render('path/to/header.hpp') == '[include_file path/to/header.hpp]'

def test_entity_declaration(entities, render):
    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.entity_declaration(entities) }}'''
    assert render(entities['cl1']) == 'class klass;'
    assert render(entities['cl3']) == textwrap.dedent('''\
        template<>
        class tmpl<int, std::string>
            : public virtual ``[link ns1__ns2__klass.nested [^ns1::ns2::klass::nested]]``;''')
    assert render(entities['e1']) == 'enum class e1\n    : unsigned;'
    assert render(entities['t1']) == 'using TA = std::string;'
    assert render(entities['v1']) == 'static\nint Var = 11;'
    assert render(entities['v2']) == 'static constexpr\ndouble Var2;'
    assert render(entities['v3']) == 'int(*Var3)(int, long);'

def test_function_declaration(entities, render):
    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.function_declaration(entities) }}'''
    assert render(entities['f1']) == textwrap.dedent('''\
        void
        func(
            int,
            int (&n) [1],
            ``[link ns1__ns2__klass [^ns1::ns2::klass]]`` = ``[link ns1__Var [^ns1::Var]]``,
            unsigned t = 2);''')
    assert render(entities['f2']) == 'static constexpr void\nfunc() noexcept(noexcept((void)12));'
    assert render(entities['cl1_c']) == 'explicit constexpr\nklass();'
    assert render(entities['cl1_d']) == '~klass() noexcept(false);'
    assert render(entities['o=']) == 'void\noperator=() noexcept;'
    assert render(entities['g1']) == 'void\ng() && noexcept(foobar(")", "("));'
    assert render(entities['g2']) == 'void\ng() const&;'
    assert render(entities['g3']) == 'void\ng() & = delete;'

def test_write_entity(cfg, entities, render):
    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.write_entity(entities) }}'''
    assert render(entities['ns1']) == textwrap.dedent('''
        [section:ns1__TA ns1::TA]
        [indexterm1 TA]


        [heading Synopsis]
        Defined in header [include_file /path/to/ns1.hpp]

        ```
        using TA = std::string;
        ```





        [endsect]


        [section:ns1__e1 ns1::e1]
        [indexterm1 e1]


        [heading Synopsis]
        Defined in header [include_file /path/to/ns1.hpp]

        ```
        enum class e1
            : unsigned;
        ```
        [heading Values]
        [table [[Name][Description]]
          [
            [[^a]
            ]
            [
            ]
          ]
        ]





        [endsect]

        [section:ns1__Var ns1::Var]
        [indexterm1 Var]


        [heading Synopsis]
        Defined in header [include_file /path/to/ns1.hpp]

        ```
        static
        int Var = 11;
        ```





        [endsect]
        [section:ns1__Var2 ns1::Var2]
        [indexterm1 Var2]


        [heading Synopsis]
        Defined in header [include_file /path/to/ns1.hpp]

        ```
        static constexpr
        double Var2;
        ```





        [endsect]
        [section:ns1__Var3 ns1::Var3]
        [indexterm1 Var3]


        [heading Synopsis]
        Defined in header [include_file /path/to/ns1.hpp]

        ```
        int(*Var3)(int, long);
        ```





        [endsect]
        ''')

    entities['cl1'].members = dict([
        (k, v) for (k, v) in entities['cl1'].members.items()
        if not v.name.startswith('o') or v.name == 'operator[]'
    ])
    cfg['default_namespace'] = 'ns1'
    cfg['convenience_header'] = 'some/header.hpp'
    assert render(entities['cl1']) == textwrap.dedent('''
        [section:ns1__ns2__klass ns2::klass]


        [heading Synopsis]
        Defined in header [include_file /path/to/ns1.hpp]

        ```
        class klass;
        ```
        [heading Types]
        [table [[Name][Description]]
          [
            [[*[link ns1__ns2__klass.nested nested]]
            ]
            [
            ]
          ]
        ]
        [heading Member Functions]
        [table [[Name][Description]]
          [
            [[*[link ns1__ns2__klass.g g]]
            ]
            [One overload of g \'''<sbr/>\'''[role silver \\u2014]\'''<sbr/>\'''Another overload of g \'''<sbr/>\'''[role silver \\u2014]\'''<sbr/>\'''Third overload of g
            ]
          ]
          [
            [[*[link ns1__ns2__klass.klass klass]\\u00A0[role silver \\[constructor\\]]]
            ]
            [
            ]
          ]
          [
            [[*[link ns1__ns2__klass.operator__lb__rb_ operator\\[\\]]]
            ]
            [This is ['one] of the operators
            ]
          ]
          [
            [[*[link ns1__ns2__klass._dtor_klass ~klass]\\u00A0[role silver \\[destructor\\]]]
            ]
            [
            ]
          ]
        ]
        [heading Static Member Functions]
        [table [[Name][Description]]
          [
            [[*[link ns1__ns2__klass.func func]]
            ]
            [
            ]
          ]
        ]
        [heading Friends]
        [table [[Name][Description]]
          [
            [[*[link ns1__ns2__klass.func_fr func]]
            ]
            [
            ]
          ]
        ]

        [heading Description]
        This class contains [link ns1__ns2__klass.nested `nested`].



        Convenience header [include_file some/header.hpp]

        [section:nested ns2::klass::nested]
        [indexterm2 nested..klass]


        [heading Synopsis]


        ```
        template<
            class,
            int (&N) [1],
            class = std::String,
            class T = std::String>
        class nested;
        ```
        [heading Types]
        [table [[Name][Description]]
          [
            [[*[link ns1__ns2__klass.nested.nested2 nested2]]
            ]
            [
            ]
          ]
        ]





        [section:nested2 ns2::klass::nested::nested2]
        [indexterm2 nested2..nested]


        [heading Synopsis]


        ```
        class nested2;
        ```

        [heading Description]
        See [link ns1__ns2__klass.nested.nested2 `nested2`].





        [endsect]



        [endsect]



        [section:g ns2::klass::g]
        [indexterm2 g..klass]
        One overload of g
        ```
        void
        ``[link ns1__ns2__klass.g.overload1 g]``() && noexcept(foobar(")", "("));
          ``[''\''&raquo;\''' [link ns1__ns2__klass.g.overload1 `more...`]]``
        ```
        Another overload of g
        ```
        void
        ``[link ns1__ns2__klass.g.overload2 g]``() const&;
          ``[''\''&raquo;\''' [link ns1__ns2__klass.g.overload2 `more...`]]``
        ```
        Third overload of g
        ```
        void
        ``[link ns1__ns2__klass.g.overload3 g]``() & = delete;
          ``[''\''&raquo;\''' [link ns1__ns2__klass.g.overload3 `more...`]]``
        ```

        [section:overload1 ns2::klass::g (1 of 3 overloads)]
        One overload of g



        [heading Synopsis]


        ```
        void
        g() && noexcept(foobar(")", "("));
        ```





        [endsect]
        [section:overload2 ns2::klass::g (2 of 3 overloads)]
        Another overload of g



        [heading Synopsis]


        ```
        void
        g() const&;
        ```





        [endsect]
        [section:overload3 ns2::klass::g (3 of 3 overloads)]
        Third overload of g



        [heading Synopsis]


        ```
        void
        g() & = delete;
        ```





        [endsect]
        [endsect]


        [section:klass ns2::klass::klass]
        [indexterm2 klass..klass]


        [heading Synopsis]


        ```
        explicit constexpr
        klass();
        ```





        [endsect]


        [section:_dtor_klass ns2::klass::~klass]
        [indexterm2 ~klass..klass]


        [heading Synopsis]


        ```
        ~klass() noexcept(false);
        ```





        [endsect]


        [section:operator__lb__rb_ ns2::klass::operator\\[\\]]
        [indexterm2 operator\\[\\]..klass]
        This is ['one] of the operators



        [heading Synopsis]


        ```
        void
        operator[]() noexcept;
        ```





        [endsect]


        [section:func ns2::klass::func]
        [indexterm2 func..klass]


        [heading Synopsis]


        ```
        static constexpr void
        func() noexcept(noexcept((void)12));
        ```





        [endsect]


        [section:func_fr ns2::klass::func]
        [indexterm2 func..klass]


        [heading Synopsis]
        Defined in header [include_file /path/to/ns1.hpp]

        ```
        void
        func(
            int,
            int (&n) [1],
            ``[link ns1__ns2__klass [^ns2::klass]]`` = ``[link ns1__Var [^Var]]``,
            unsigned t = 2);
        ```



        Convenience header [include_file some/header.hpp]

        [endsect]



        [endsect]

        ''')

def test_defaulted(cfg, render):
    func = {
        'tag': 'memberdef',
        'kind': 'function',
        'id': 'func',
        'items': [
            { 'tag': 'name', 'items': ['func'] },
            { 'tag': 'argsstring', 'items': ['() =default'] },
            { 'tag': 'type', 'items': ['void'] },
        ]
    }
    ns = docca.Namespace(
        make_elem({
            'tag': 'compound',
            'id': 'ns',
            'items': [
                { 'tag': 'compoundname', 'items': ['ns'] },
                { 'tag': 'sectiondef', 'items': [func] },
            ]
        }))
    ns.resolve_references()
    func = ns.index['func']
    func.resolve_references()

    render.template = '''\
        {%- import "docca/quickbook/components.jinja2" as qbk -%}
        {{ qbk.function_declaration(entities) }}'''
    assert render(func) == textwrap.dedent('''\
        void
        func();''')

    render.env.globals['Config']['show_defaulted'] = True
    assert render(func) == textwrap.dedent('''\
        void
        func() = default;''')
