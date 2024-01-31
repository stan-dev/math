#!/usr/bin/python3
#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# Generates the seed corpus from captured wireshark frames.

import json
from enum import Enum
from pathlib import Path
from os import path
import struct
from typing import List, NamedTuple
import tarfile
from io import BytesIO

_REPO_BASE = Path(path.join(path.dirname(path.realpath(__file__)), '..', '..')).absolute()
_TEXT_JSON = _REPO_BASE.joinpath('private', 'fuzzing', 'text.json')
_BINARY_JSON = _REPO_BASE.joinpath('private', 'fuzzing', 'binary.json')

class _ProtocolFieldType(Enum):
    decimal = 0x00
    tiny = 0x01
    short_ = 0x02
    long_ = 0x03
    float_ = 0x04
    double_ = 0x05
    null = 0x06
    timestamp = 0x07
    longlong = 0x08
    int24 = 0x09
    date = 0x0a
    time = 0x0b
    datetime = 0x0c
    year = 0x0d
    varchar = 0x0f
    bit = 0x10
    json = 0xf5
    newdecimal = 0xf6
    blob = 0xfc
    var_string = 0xfd
    string = 0xfe
    geometry = 0xff


class _ColumnType(Enum):
    tinyint = 0
    smallint = 1
    mediumint = 2
    int_ = 3
    bigint = 4
    float_ = 5
    double_ = 6
    decimal = 7
    bit = 8
    year = 9
    time = 10
    date = 11
    datetime = 12
    timestamp = 13
    char_ = 14
    varchar = 15
    binary = 16
    varbinary = 17
    text = 18
    blob = 19
    enum_ = 20
    set = 21
    json = 22
    geometry = 23
    unknown = 24


class _Encoding(Enum):
    text = 0
    binary = 1


_TYPE_MAP = {
    _ProtocolFieldType.decimal: _ColumnType.decimal,
    _ProtocolFieldType.newdecimal: _ColumnType.decimal,
    _ProtocolFieldType.geometry: _ColumnType.geometry,
    _ProtocolFieldType.tiny: _ColumnType.tinyint,
    _ProtocolFieldType.short_: _ColumnType.smallint,
    _ProtocolFieldType.int24: _ColumnType.mediumint,
    _ProtocolFieldType.long_: _ColumnType.int_,
    _ProtocolFieldType.longlong: _ColumnType.bigint,
    _ProtocolFieldType.float_: _ColumnType.float_,
    _ProtocolFieldType.double_: _ColumnType.double_,
    _ProtocolFieldType.bit: _ColumnType.bit,
    _ProtocolFieldType.date: _ColumnType.date,
    _ProtocolFieldType.datetime: _ColumnType.datetime,
    _ProtocolFieldType.timestamp: _ColumnType.timestamp,
    _ProtocolFieldType.time: _ColumnType.time,
    _ProtocolFieldType.year: _ColumnType.year,
    _ProtocolFieldType.json: _ColumnType.json,
}


def _compute_column_type(mysql_frame) -> _ColumnType:
    type = _ProtocolFieldType(int(mysql_frame['mysql.field.type_raw'][0], 16))
    isbin = mysql_frame['mysql.field.charsetnr_raw'][0] == '3f00'
    if type == _ProtocolFieldType.string:
        if mysql_frame['mysql.field.flags_tree']['mysql.field.flags.set_raw'][0] != '0':
            return _ColumnType.set
        elif mysql_frame['mysql.field.flags_tree']['mysql.field.flags.enum_raw'][0] != '0':
            return _ColumnType.enum_
        elif isbin:
            return _ColumnType.binary
        else:
            return _ColumnType.char_
    elif type == _ProtocolFieldType.var_string:
        if isbin:
            return _ColumnType.varbinary
        else:
            return _ColumnType.varchar
    elif type == _ProtocolFieldType.blob:
        if isbin:
            return _ColumnType.blob
        else:
            return _ColumnType.text
    else:
        return _TYPE_MAP[type]


class _Sample(NamedTuple):
    name: str
    content: bytes


def _encoding_to_str(enc: _Encoding) -> str:
    return 'text' if enc == _Encoding.text else 'binary'


# meta: spans 2 bytes
# meta[0][low 7 bits]: column_type
# meta[0][high bit]: is unsigned flag
# meta[1]: decimals
def _gen_meta(mysql_frame) -> bytes:
    type = _compute_column_type(mysql_frame)
    unsigned_flag = int(mysql_frame['mysql.field.flags_tree']['mysql.field.flags.unsigned_raw'][0])
    assert unsigned_flag == 0 or unsigned_flag == 1
    decimals = int(mysql_frame['mysql.field.decimals_raw'][0], 16)
    return struct.pack(
        '<BB',
        type.value | (unsigned_flag << 7),
        decimals
    )

# Parses a text row into its fields. Relies on int_lenencs not being more than 1 byte
def _parse_text_row(r: bytes, num_fields: int) -> List[bytes]:
    res: List[bytes] = []
    offset = 0
    for _ in range(num_fields):
        (length,) = struct.unpack('<B', r[offset:offset+1])
        if length >= 0xfb: # Too long or NULL samples
            continue
        res.append(r[offset+1:offset+1+length])
        offset += length + 1
    return res


def _gen_text_field_samples(response, table: str) -> List[_Sample]:
    mysql_raw = response['_source']['layers']['mysql_raw']
    mysql = response['_source']['layers']['mysql']

    num_fields = int(mysql[0]['mysql.num_fields_raw'][0], 16)
    metas = [_gen_meta(mysql[i + 1]) for i in range(num_fields)]

    res: List[_Sample] = []
    for i, rec in enumerate(mysql_raw[1 + num_fields:-1]):
        row = bytes.fromhex(rec[0][2*4:])
        fields = _parse_text_row(row, num_fields)
        res += [_Sample(
            f'fuzz_text_field/{table}_{i}_{j}.bin',
            meta + field
        ) for j, (meta, field) in enumerate(zip(metas, fields))]

    return res


def _gen_row_samples(response, enc: _Encoding, table: str) -> List[_Sample]:
    mysql_raw = response['_source']['layers']['mysql_raw']
    mysql = response['_source']['layers']['mysql']

    # Header[0][low 7 bits]: num_fields
    # Header[0][high bit]: encoding
    num_fields = int(mysql[0]['mysql.num_fields_raw'][0], 16)
    header = struct.pack('<B', num_fields | (enc.value << 7))

    # As many meta blocks as num_fields
    meta = b''.join(_gen_meta(mysql[i + 1]) for i in range(num_fields))

    # Rows
    return [
        _Sample(
            f'fuzz_row/{_encoding_to_str(enc)}_{table}_{i}.bin',
            b''.join((header, meta, bytes.fromhex(rec[0][2*4:])))
        )
        for i, rec in enumerate(mysql_raw[1 + num_fields:-1])
    ]


# Frames should be captured with wireshark (see below for each protocol encoding)
# and saved as pcapng files. They should be converted to jsonraw using
# tshark -r fname.pcapng -T jsonraw --no-duplicate-keys -Y 'mysql' > fname.json
def main():
    corpus: List[_Sample] = []

    # Text. These should contain frames generated by running SELECT * FROM <table>
    # repeatedly over the set of tables, from mysql command line.
    # This generates two frames (request and response) for each select
    with open(_TEXT_JSON, 'rt') as f:
        obj = json.load(f)
    for query, response in zip(*(iter(obj),) * 2):
        table = bytes.fromhex(query['_source']['layers']['mysql']['mysql.payload_raw'][0])[1:].decode().split(' ')[-1]
        corpus += _gen_row_samples(response, _Encoding.text, table)
        corpus += _gen_text_field_samples(response, table)

    # Binary. These should contain frames generated by running SELECT * FROM <table>
    # repeatedly over the set of tables, from the Python mysql protocol driver,
    # with code like:
    #    cursor = cnx.cursor(prepared=True)
    #    cursor.execute(f'SELECT * FROM {table}')
    #    cursor.fetchall()
    # This generates 7 frames for each table, since the driver prepares, executes and resets the statement.
    with open(_BINARY_JSON, 'rt') as f:
        obj = json.load(f)
    for reqclose, query, qres, reset, resetres, exec, response in zip(*(iter(obj),) * 7):
        table = bytes.fromhex(query['_source']['layers']['mysql']['mysql.request_raw'][0])[1:].decode().split(' ')[-1]
        corpus += _gen_row_samples(response, _Encoding.binary, table)
    
    # Frames for simpler deserialization functions are stored in a json file
    with open(_REPO_BASE.joinpath('tools', 'scripts', 'seed_corpus_samples.json'), 'rt') as f:
        obj = json.load(f)
    for fuzzer, samples in obj.items():
        corpus += [_Sample(
            f'{fuzzer}/{sample["name"]}.bin',
            bytes.fromhex(sample["msg"])
        ) for sample in samples]
    
    # Write the final file
    corpus_path = _REPO_BASE.joinpath('test', 'fuzzing', 'seedcorpus.tar.gz')
    with tarfile.open(corpus_path, "w:gz") as tar:
        for sample in corpus:
            tinfo = tarfile.TarInfo(sample.name)
            tinfo.size = len(sample.content)
            tar.addfile(tinfo, BytesIO(sample.content))


if __name__ == '__main__':
    main()
