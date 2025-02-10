#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import os
from pathlib import Path
from os import path
from typing import List, Tuple
import glob
from abc import abstractmethod, ABCMeta
import sys

THIS_FOLDER = path.abspath(path.dirname(path.realpath(__file__)))
sys.path.append(path.join(THIS_FOLDER))

import examples_qbk

# Script to get file headers (copyright notices
# and include guards) okay and up to date

VERBOSE = False

REPO_BASE = path.abspath(path.join(path.dirname(path.realpath(__file__)), '..', '..'))
BASE_FOLDERS = [
    'cmake',
    'doc',
    'example',
    'include',
    'src',
    'test',
    'tools',
    'bench',
    '.github'
]
BASE_FILES = [
    'CMakeLists.txt',
    '.drone.star',
    '.codecov.yml'
]
HTML_GEN_PATH = path.join(REPO_BASE, 'doc', 'html')

HEADER_TEMPLATE = '''{begin}
{linesym} Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
{linesym}
{linesym} Distributed under the Boost Software License, Version 1.0. (See accompanying
{linesym} file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
{end}'''

SRC_HPP_TEMPLATE = '''
// This file is meant to be included once, in a translation unit of
// the program, with the macro BOOST_MYSQL_SEPARATE_COMPILATION defined.

#include <boost/mysql/detail/config.hpp>

#ifndef BOOST_MYSQL_SEPARATE_COMPILATION
#error You need to define BOOST_MYSQL_SEPARATE_COMPILATION in all translation units that use the compiled version of Boost.MySQL, \\
    as well as the one where this file is included.
#endif

{includes}

#endif'''

MYSQL_ERROR_HEADER = path.join(REPO_BASE, 'private', 'mysqld_error.h')
MARIADB_ERROR_HEADER = path.join(REPO_BASE, 'private', 'mariadb_error.h')

def find_first_blank(lines: List[str]) -> int:
    return next((i for i, line in enumerate(lines) if line.replace('\n', '') == ''), len(lines))

def read_file(fpath):
    with open(fpath, 'rt') as f:
        try:
            return f.readlines()
        except Exception as err:
            raise SystemError(f'Error processing file {fpath}') from err
    
def write_file(fpath, lines):
    with open(fpath, 'wt') as f:
        f.writelines(lines)

def text_to_lines(text):
    return [line + '\n' for line in text.split('\n')]

def gen_header(linesym, opensym=None, closesym=None, shebang=None, include_guard=None):
    opensym = linesym if opensym is None else opensym
    closesym = linesym if closesym is None else closesym
    if shebang is None:
        begin = opensym
    else:
        begin = shebang + '\n' + opensym
    if include_guard is None:
        end = closesym
    else:
        end = closesym + '\n\n#ifndef {0}\n#define {0}'.format(include_guard)
    return text_to_lines(HEADER_TEMPLATE.format(begin=begin, end=end, linesym=linesym))

class BaseProcessor(metaclass=ABCMeta):
    skip = False

    @abstractmethod
    def process(self, lines: List[str], fpath: str) -> List[str]:
        return lines
    
    name = ''

class NormalProcessor(BaseProcessor):
    def __init__(self, name, header):
        self.header = header
        self.name = name
        
    def process(self, lines: List[str], _: str) -> List[str]:
        first_blank = find_first_blank(lines)
        lines = self.header + lines[first_blank:]
        return lines
        
class HppProcessor(BaseProcessor):
    name = 'hpp'
    
    def process(self, lines: List[str], fpath: str) -> List[str]:
        first_content = [i for i, line in enumerate(lines) if line.startswith('#define')][0] + 1
        iguard = self._gen_include_guard(fpath)
        header = gen_header('//', include_guard=iguard)
        lines = header + lines[first_content:]
        return lines
        
        
    @staticmethod
    def _gen_include_guard(fpath):
        include_base = path.join(REPO_BASE, 'include')
        if fpath.startswith(include_base):
            relpath = path.relpath(fpath, include_base)
        else:
            relpath = path.join('boost', 'mysql', path.relpath(fpath, REPO_BASE))
        return relpath.replace('/', '_').replace('.', '_').upper()


class SrcHppProcessor(HppProcessor):
    name = 'src.hpp'

    def process(self, _: List[str], fpath: str) -> List[str]:
        base = Path(REPO_BASE)
        includes = [
            fname.relative_to(base.joinpath('include'))
            for fname in sorted(base.joinpath('include', 'boost', 'mysql', 'impl').rglob('*.ipp'))
        ]
        return gen_header('//', include_guard=self._gen_include_guard(fpath)) + \
            text_to_lines(
                SRC_HPP_TEMPLATE.format(
                    includes='\n'.join(f'#include <{inc}>' for inc in includes)
                )
            )


class MysqlHppProcessor(HppProcessor):
    name = 'mysql.hpp'

    def process(self, _: List[str], fpath: str) -> List[str]:
        base = Path(REPO_BASE)
        exclusions = ('src.hpp', 'pfr.hpp')
        includes = [
            fname.relative_to(base.joinpath('include'))
            for fname in sorted(base.joinpath('include', 'boost', 'mysql').glob('*.hpp'))
            if fname.name not in exclusions
        ]
        return gen_header('//', include_guard=self._gen_include_guard(fpath)) + \
            ['\n'] + \
            [f'#include <{inc}>\n' for inc in includes] + \
            ['\n', '#endif\n']


class XmlProcessor(BaseProcessor):
    name = 'xml'
    header = gen_header('   ', '<!--', '-->')
    
    def process(self, lines: List[str], fpath: str) -> List[str]:
        if lines[0].startswith('<?'):
            first_blank = [i for i, line in enumerate(lines) if line.strip() == ''][0]
            first_content = [i for i, line in enumerate(lines[first_blank:]) \
                             if line.startswith('<') and not line.startswith('<!--')][0] + first_blank
            lines = lines[0:first_blank] + ['\n'] + self.header + ['\n'] + lines[first_content:]
        else:
            lines = NormalProcessor('xml', self.header).process(lines, fpath)
        
        return lines
        
        
class IgnoreProcessor(BaseProcessor):
    name = 'ignore'
    skip = True
    
    def process(self, lines: List[str], _: str) -> List[str]:
        return lines
        
hash_processor = NormalProcessor('hash', gen_header('#'))
qbk_processor = NormalProcessor('qbk', gen_header('   ', opensym='[/', closesym=']'))
sql_processor = NormalProcessor('sql', gen_header('--'))
cpp_processor = NormalProcessor('cpp', gen_header('//'))
py_processor = NormalProcessor('py', gen_header('#', shebang='#!/usr/bin/python3'))
bash_processor = NormalProcessor('bash', gen_header('#', shebang='#!/bin/bash'))
bat_processor = NormalProcessor('bat', gen_header('@REM'))

FILE_PROCESSORS : List[Tuple[str, BaseProcessor]] = [
    ('src.hpp', SrcHppProcessor()),
    ('mysql.hpp', MysqlHppProcessor()),
    ('CMakeLists.txt', hash_processor),
    ('.cmake', hash_processor),
    ('.cmake.in', hash_processor),
    ('Jamfile', hash_processor),
    ('.jam', hash_processor),
    ('Doxyfile', hash_processor),
    ('.qbk', qbk_processor),
    ('.sql', sql_processor),
    ('.py', py_processor),
    ('.sh', bash_processor),
    ('.bat', bat_processor),
    ('.ps1', hash_processor),
    ('.yml', hash_processor),
    ('.cnf', hash_processor),
    ('.dockerfile', hash_processor),
    ('.star', hash_processor),
    ('.cpp', cpp_processor),
    ('.hpp', HppProcessor()),
    ('.ipp', HppProcessor()),
    ('.xml', XmlProcessor()),
    ('.xsl', XmlProcessor()),
    ('.svg', IgnoreProcessor()),
    ('valgrind_suppressions.txt', IgnoreProcessor()),
    ('.pem', IgnoreProcessor()),
    ('.md', IgnoreProcessor()),
    ('.csv', IgnoreProcessor()),
    ('.tar.gz', IgnoreProcessor()),
    ('.json', IgnoreProcessor()),
    ('.txt', IgnoreProcessor()),
    ('.pyc', IgnoreProcessor()),
]

def process_file(fpath: str):
    try:
        for ext, processor in FILE_PROCESSORS:
            if fpath.endswith(ext):
                if VERBOSE:
                    print('Processing file {} with processor {}'.format(fpath, processor.name))
                if not processor.skip:
                    lines = read_file(fpath)
                    output_lines = processor.process(lines, fpath)
                    if output_lines != lines:
                        write_file(fpath, output_lines)
                break
        else:
            raise ValueError('Could not find a suitable processor for file: ' + fpath)
    except Exception:
        print(f'Found error processing {fpath}')
        raise
    
def process_all_files():
    for base_folder in BASE_FOLDERS:
        base_folder_abs = path.join(REPO_BASE, base_folder)
        for curdir, _, files in os.walk(base_folder_abs):
            if curdir.startswith(HTML_GEN_PATH):
                if VERBOSE:
                    print('Ignored directory {}'.format(curdir))
                continue
            for fname in files:
                process_file(path.join(curdir, fname))
    for fname in BASE_FILES:
        process_file(path.join(REPO_BASE, fname))


# Check that cmake and b2 test source files are equal
def verify_test_consistency():
    for test_type in ('unit', 'integration'):
        for ftocheck in ('Jamfile', 'CMakeLists.txt'):
            base_path = path.join(REPO_BASE, 'test', test_type)
            tests = glob.glob(base_path + '/**/*.cpp', recursive=True)
            tests = [elm.replace(base_path + '/', '') for elm in tests]

            with open(path.join(REPO_BASE, 'test', test_type, ftocheck), 'rt') as f:
                contents = f.readlines()
                contents = ''.join(elm for elm in contents if not '#' in elm)

            for t in tests:
                if ' ' + t not in contents:
                    print(f'File {t} not in {test_type}/{ftocheck}')

            
def main():
    process_all_files()
    verify_test_consistency()
    examples_qbk.main()
            
            
if __name__ == '__main__':
    main()
        
