#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from pathlib import Path
from os import makedirs
from io import BytesIO
import struct
from datetime import datetime, date
from typing import NamedTuple, List
import csv
import itertools
import shutil
from .common import REPO_BASE


SEED_CORPUS_PATH = Path('/tmp/seedcorpus')


class _Sample(NamedTuple):
    fuzzer: str
    name: str
    content: bytes


# Generates a list of samples by reading one of the CSV files
def _read_csv(fname: str) -> List[_Sample]:
    csv_path = REPO_BASE.joinpath('tools', 'seed_corpus', fname)
    with open(csv_path, 'rt') as f:
        return [
            _Sample(s['fuzzer'], s['name'], bytes.fromhex(s['content']))
            for s in csv.DictReader(f)
        ]


def _gen_utf8mb4_valid() -> List[_Sample]:
    samples = [
        ("euro_symbol", "€."),
        ("greek", "Μπορώ να φάω σπασμένα γυαλιά χωρίς να πάθω τίποτα."),
        ("icelandic", "Ég get etið gler án þess að meiða mig."),
        ("polish", "Mogę jeść szkło, i mi nie szkodzi."),
        ("romanian", "Pot să mănânc sticlă și ea nu mă rănește."),
        ("ukrainian", "Я можу їсти шкло, й воно мені не пошкодить."),
        ("armenian", "Կրնամ ապակի ուտել և ինծի անհանգիստ չըներ։"),
        ("georgian", "მინას ვჭამ და არა მტკივა."),
        ("hindi", "मैं काँच खा सकता हूँ, मुझे उस से कोई पीडा नहीं होती."),
        ("hebrew", "אני יכול לאכול זכוכית וזה לא מזיק לי."),
        ("yiddish", "איך קען עסן גלאָז און עס טוט מיר נישט װײ."),
        ("arabic", "أنا قادر على أكل الزجاج و هذا لا يؤلمني."),
        ("japanese", "私はガラスを食べられます。それは私を傷つけません。"),
        ("thai", "ฉันกินกระจกได้ แต่มันไม่ทำให้ฉันเจ็บ"),
        ("emoji", "😋😛😜🤪"),
    ]

    return [
        _Sample('fuzz_utf8mb4', name, content.encode())
        for name, content in samples
    ]


def _gen_utf8mb4_invalid() -> List[_Sample]:
    samples = [
        ("c2_ltmin", "c27f"),
        ("c2_gtmax", "c2c0"),
        ("de_gtmax", "dec0"),
        ("df_ltmin", "df7f"),
        ("df_gtmax", "dfc0"),
        ("e0_ltmin_ok", "e09f91"),
        ("e3_ltmin_ok", "e37f91"),
        ("ed_ok_ltmin", "eda07f"),
        ("ed_surrogate_min", "eda080"),
        ("f0_ltmin_ok_ok", "f08f8080"),
        ("f2_ok_gtmax_ok", "f2a1c0a3"),
        ("f4_ok_ok_gtmax", "f4a1a2c0"),
        ("overlong_slash_2byte", "c0af"),
        ("overlong_slash_3byte", "e080af"),
        ("overlong_slash_4byte", "f08080af"),
        ("overlong_slash_5byte", "f8808080af"),
        ("overlong_slash_6byte", "f880808080af"),
    ]

    return [
        _Sample('fuzz_utf8mb4', name, bytes.fromhex(content))
        for name, content in samples
    ]


def _gen_escape_string() -> List[_Sample]:
    # Some base cases that cover most of the escaping we do
    base_cases = [
        ("empty",                  b""),
        ("no_escape_ascii",        b"this is A test string"),
        ("escape_quotes",          b'some "dq"\'sq\' `bt` \\bl\\'),
        ("escape_others",          b"this '\0' null \x1a \n\r some chrs \\"),
        ("utf8_2byte",              "2byte \" ñ UTF-8\\ ò` \\".encode()),
        ("utf8_3byte",              "3byte '\uffff UTF-8'".encode()),
        ("utf8_4byte",              "4byte \r'𐀀 UTF-8\n".encode()),
        ("invalid_utf8",           b"This \"has\" invalid \xc3\\ chars"),
        ("injection_1",            b"' or \""),
        ("injection_2",            b"' OR 'x'='x\\"),
        ("injection_3",            b"'''''''''''''UNION SELECT '2"),
    ]

    # The first byte in the input contains options for backslash_escapes
    # and the quoting context
    backslash_escapes = [
        (0,      'bn'),
        (1 << 0, 'by'),
    ]

    quot_ctx = [
        (0,               'squot'),
        (1 << 1,          'btick'),
        (1 << 1 | 1 << 2, 'dquot'),
    ]

    # Generate all the cases by dot product
    cases: List[_Sample] = []
    for bl_byte, bl_name in backslash_escapes:
        for quot_byte, quot_name in quot_ctx:
            cases += [_Sample(
                'fuzz_escape_string',
                f'{bl_name}_{quot_name}_{name}',
                struct.pack('<B', bl_byte | quot_byte) + value
            ) for name, value in base_cases]
        
    return cases


def _gen_format_string() -> List[_Sample]:
    cases = [
        ("no_replacement", "1"),
        ("escaped_curly_1", "'{{}}'"),
        ("escaped_curly_2", "'}}}}{{'"),
        ("escaped_curly_3", "'}}'"),
        ("escaped_curly_4", "'{{name}}'"),
        ("one_replacement", "{} OR 1=1"),
        ("replacement_end", "{}"),
        ("replacement_start", "{} AND {}"),
        ("several_replacements", "{}, {}, {}, {}, {}, {}, {}, {}"),
        ("indexed_args_1", "{0}"),
        ("indexed_args_2", "{1}, {0}"),
        ("indexed_args_3", "{07}, {1}, {0}, {2}"),
        ("named_args_1", "{name}, {val}"),
        ("named_args_2", "{_name}, {0}, {k}, {v}"),
        ("named_args_3", "{}, {name}, {}, {k}"),
        ("utf8_1", "`eñe` + {};"),
        ("utf8_2", "uni 🤡 {} code"),
        ("unbalanced_1", "{} { bad"),
        ("unbalanced_2", "}}} bad"),
        ("invalid_name_1", "{0name}"),
        ("invalid_name_2", "{ name }"),
        ("invalid_name_3", "{eñe}"),
        ("invalid_index_1", "{999999}"),
        ("invalid_index_2", "{0x10}"),
        ("invalid_index_3", "{4.2}"),
        ("index_to_manual", "{0}, {}"),
        ("auto_to_manual", "{}, {1}, {2}"),
    ]

    return [_Sample(
        'fuzz_format_strings',
        name,
        content.encode()
    ) for name, content in cases]


def _gen_format_args() -> List[_Sample]:
    # Helper to output samples
    class Packer:
        def __init__(self) -> None:
            self.io = BytesIO()
        
        def add_null(self, v: None):
            pass
        
        def add_int64(self, v: int):
            self.io.write(struct.pack('<q', v))
        
        def add_uint64(self, v: int):
            self.io.write(struct.pack('<Q', v))
        
        def add_float(self, v: float):
            self.io.write(struct.pack('<f', v))
        
        def add_double(self, v: float):
            self.io.write(struct.pack('<d', v))

        def add_string(self, v: str):
            self.add_blob(v.encode())
        
        def add_blob(self, v: bytes):
            self.io.write(struct.pack('<B', len(v)))
            self.io.write(v)
        
        def add_date(self, v: date):
            self.io.write(struct.pack('<HBB', v.year, v.month, v.day))
        
        def add_datetime(self, v: datetime):
            self.io.write(struct.pack('<HBBBBBL', v.year, v.month, v.day, v.hour, v.minute, v.second, v.microsecond))
        
        def pack(self, type1: int, type2: int, fn1, fn2, val1, val2) -> bytes:
            # Sample format: type code, value 1, value 2
            self.io.write(struct.pack('<B', type1 | type2 << 4))
            fn1(self, val1)
            fn2(self, val2)
            return self.io.getvalue()

    # Possible types, with type codes
    types = [
        ('null',     0x00, Packer.add_null, None, None),
        ('int64',    0x01, Packer.add_int64, -1, 42),
        ('uint64',   0x02, Packer.add_uint64, 0xffffffffff, 23),
        ('float',    0x03, Packer.add_float, 4.2, -10e20),
        ('double',   0x04, Packer.add_double, -2.1e-216, 0.0),
        ('string',   0x05, Packer.add_string, 'ab\0\n\'\\"', 'm\r`\\\\abc'),
        ('blob',     0x06, Packer.add_blob,   b'\0ab\\\n`', b'a'*64),
        ('date',     0x07, Packer.add_date, date(2021, 10, 11), date(1970, 1, 1)),
        ('datetime', 0x08, Packer.add_datetime, datetime(2021, 11, 9, 10, 1, 20, 91), datetime(2100, 10, 1)),
        ('time',     0x09, Packer.add_int64, -90, 439389289202),
    ]

    # Perform a dot product of all cases
    cases: List[_Sample] = []
    for name1, type1, fn1, val1, _ in types:
        cases += [
            _Sample(
                'fuzz_format_args',
                f'{name1}_{name2}',
                Packer().pack(type1, type2, fn1, fn2, val1, val2)
            )
            for name2, type2, fn2, _, val2 in types
        ]
    
    return cases

def _gen_format_identifiers() -> List[_Sample]:
    cases = [
        '',
        'employee',
        'some_table',
        'with spaces',
        'wíth uñicode',
        'need`ds esc``aping',
        '`',
        '\'"ab[]def',
    ]

    return [
        _Sample('fuzz_format_identifier', str(i), s.encode())
        for i, s in enumerate(cases)
    ]


# The payloads in the text file have been generated by running
# sqlmap (https://sqlmap.org/) on a test server
def _gen_format_sql_injection() -> List[_Sample]:
    with open(REPO_BASE.joinpath('tools', 'seed_corpus', 'sql_injection_payloads.txt'), 'rt', encoding='utf-8') as f:
        payloads = filter(lambda p: p != '', f.read().split('\n'))
    
    return [_Sample(
        'fuzz_format_sql_injection',
        str(i),
        payload.encode()
    ) for i, payload in enumerate(payloads)]


def generate_seed_corpus():
    # Generate the samples in memory
    samples = list(itertools.chain(
        _read_csv('field_table.csv'),
        _read_csv('protocol_messages.csv'),
        _gen_utf8mb4_valid(),
        _gen_utf8mb4_invalid(),
        _gen_escape_string(),
        _gen_format_string(),
        _gen_format_args(),
        _gen_format_sql_injection(),
        _gen_format_identifiers(),
    ))

    # Generate the directory structure
    fuzzers = list(set(s.fuzzer for s in samples))
    if SEED_CORPUS_PATH.exists():
        shutil.rmtree(SEED_CORPUS_PATH, )
    for f in fuzzers:
        makedirs(SEED_CORPUS_PATH.joinpath(f))
    
    # Generate the samples
    for s in samples:
        with open(SEED_CORPUS_PATH.joinpath(s.fuzzer, s.name + '.bin'), 'wb') as f:
            f.write(s.content)
