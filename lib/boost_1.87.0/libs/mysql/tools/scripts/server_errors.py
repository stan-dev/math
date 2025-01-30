#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

# This scripts generates files containing the server defined error codes,
# and code to convert from error codes to strings. This is complex because:
#  - There are *a lot* of error codes.
#  - There are common error codes and MariaDB/MySQL specific ones.
#  - Some codes have been repurposed, renamed or removed from MySQL 5.x to MySQL 8.x and MariaDB.
# To generate precise output, we need the mysqld_error.h header for MySQL 5.x, 8.x and MariaDB.
# MySQL and MariaDB often rename, deprecate, etc. codes. To attempt to maintain sanity,
# we generate a CSV file with what we did for the previous version. We load it and combine it
# with the new errors, and we only perform backwards-compatible changes.
#
# To update the codes for a new release:
#   - MySQL: https://dev.mysql.com/downloads/mysql/
#   - MariaDB: https://mariadb.org/connector-c/all-releases/

import pandas as pd
from os import path
from pathlib import Path
from typing import Literal, List, Optional, cast, NamedTuple
from subprocess import run


# DataFrames have 'symbol', 'numbr' columns
class ErrorCodes(NamedTuple):
    common: pd.DataFrame
    mysql: pd.DataFrame
    mariadb: pd.DataFrame


REPO_BASE = Path(path.abspath(path.join(path.dirname(path.realpath(__file__)), '..', '..')))
_CSV_PATH = REPO_BASE.joinpath('tools', 'error_codes.csv')

# All server errors range between 1000 and 4999. Errors between 2000 and 2999
# are reserved for the client and are not used. In theory, codes between COMMON_ERROR_FIRST
# and COMMON_ERROR_LAST are shared between MySQL and MariaDB. However, some exceptions apply -
# some codes were not used originally by MySQL and are now used only by MariaDB, some have been renamed,
# etc. All codes >= COMMON_ERROR_LAST are server-specific. Codes between [COMMON_ERROR_FIRST, COMMON_ERROR_LAST)
# may be euther common or server-specific.
COMMON_ERROR_FIRST = 1000
COMMON_ERROR_LAST = 1880
SERVER_ERROR_LAST = 5000

COMMON_SERVER_ERRC_ENTRY = '''
    /**
     * \\brief Common server error. Error number: {number}, symbol:
     * <a href="https://dev.mysql.com/doc/mysql-errors/8.0/en/server-error-reference.html#error_{symbol_lower}">{symbol_upper}</a>.
     */
    {symbol_lower} = {number},
'''

COMMON_SERVER_ERRC_TEMPLATE = '''
#ifndef BOOST_MYSQL_COMMON_SERVER_ERRC_HPP
#define BOOST_MYSQL_COMMON_SERVER_ERRC_HPP

#include <boost/mysql/error_code.hpp>

#include <boost/mysql/detail/config.hpp>

#include <boost/system/error_category.hpp>

namespace boost {{
namespace mysql {{

/**
 * \\brief Server-defined error codes, shared between MySQL and MariaDB.
 * \\details The numeric value and semantics match the ones described in the MySQL documentation.
 * For more info, consult the error reference for
 * <a href="https://dev.mysql.com/doc/mysql-errors/8.0/en/server-error-reference.html">MySQL 8.0</a>, 
 * <a href="https://dev.mysql.com/doc/mysql-errors/5.7/en/server-error-reference.html">MySQL 5.7</a>,
 * <a href="https://mariadb.com/kb/en/mariadb-error-codes/">MariaDB</a>.
 */
enum class common_server_errc : int
{{
{}
}};

BOOST_MYSQL_DECL
const boost::system::error_category& get_common_server_category() noexcept;

/// Creates an \\ref error_code from a \\ref common_server_errc.
inline error_code make_error_code(common_server_errc error)
{{
    return error_code(static_cast<int>(error), get_common_server_category());
}}

}}  // namespace mysql

#ifndef BOOST_MYSQL_DOXYGEN
namespace system {{

template <>
struct is_error_code_enum<::boost::mysql::common_server_errc>
{{
    static constexpr bool value = true;
}};

}}  // namespace system
#endif

}}  // namespace boost

#ifdef BOOST_MYSQL_HEADER_ONLY
#include <boost/mysql/impl/error_categories.ipp>
#endif

#endif
'''

# Render the enumeration with common codes
def render_common_server_errc(df_common: pd.DataFrame) -> str:
    entries = ''.join(COMMON_SERVER_ERRC_ENTRY.format(
        number=r.numbr,
        symbol_upper=r.symbol,
        symbol_lower=r.symbol.lower()
    ) for r in df_common.itertuples())
    return COMMON_SERVER_ERRC_TEMPLATE.format(entries)


SPECIFIC_SERVER_ERRC_ENTRY = '''
/// Server error specific to {flavor}. Error number: {number}, symbol: {symbol_upper}.
BOOST_INLINE_CONSTEXPR int {symbol_lower} = {number};
'''

SPECIFIC_SERVER_ERRC_TEMPLATE = '''
#ifndef BOOST_MYSQL_{flavor}_SERVER_ERRC_HPP
#define BOOST_MYSQL_{flavor}_SERVER_ERRC_HPP

#include <boost/config.hpp>

namespace boost {{
namespace mysql {{

namespace {flavor}_server_errc {{
{entries}
}}  // namespace {flavor}_server_errc

}}  // namespace mysql
}}  // namespace boost

#endif
'''

# Render error codes specific to a server
def render_server_specific_errc(flavor: Literal['mysql', 'mariadb'], df_db: pd.DataFrame):
    entries = ''.join(SPECIFIC_SERVER_ERRC_ENTRY.format(
        number=r.numbr,
        symbol_upper=r.symbol,
        symbol_lower=r.symbol.lower(),
        flavor=flavor
    ) for r in df_db.itertuples())
    return SPECIFIC_SERVER_ERRC_TEMPLATE.format(entries=entries, flavor=flavor)

SERVER_ERROR_TO_STRING_TEMPLATE_ENTRY='    case {number}: return "{symbol_lower}";\n'

SERVER_ERROR_TO_STRING_TEMPLATE='''
#ifndef BOOST_MYSQL_IMPL_INTERNAL_ERROR_SERVER_ERROR_TO_STRING_HPP
#define BOOST_MYSQL_IMPL_INTERNAL_ERROR_SERVER_ERROR_TO_STRING_HPP

// This file was generated by server_errors.py - do not edit directly.

#include <boost/config.hpp>
#include <boost/mysql/impl/internal/error/server_error_to_string.hpp>

namespace boost {{
namespace mysql {{
namespace detail {{

BOOST_INLINE_CONSTEXPR const char* common_error_messages[] = {{
{common_entries}
}};

}} // namespace detail
}} // namespace mysql
}} // namespace boost

const char* boost::mysql::detail::common_error_to_string(int v)
{{
    constexpr int first = {common_error_first};
    constexpr int last = first + sizeof(common_error_messages) / sizeof(const char*);
    return (v >= first && v < last) ? common_error_messages[v - first] : nullptr;
}}

const char* boost::mysql::detail::mysql_error_to_string(int v)
{{
    switch (v)
    {{
{mysql_entries}
    default: return "<unknown MySQL-specific server error>";
    }}
}}

const char* boost::mysql::detail::mariadb_error_to_string(int v)
{{
    switch (v)
    {{
{mariadb_entries}
    default: return "<unknown MariaDB-specific server error>";
    }}
}}

#endif

'''

# Renders the cpp that implements server error codes to strings
def render_server_error_to_string(df_common: pd.DataFrame, df_mysql: pd.DataFrame, df_mariadb: pd.DataFrame) -> str:
    # Common entries. We need to include non-present entries here, too (as nullptr's)
    number_to_symbol = df_common.set_index('numbr')['symbol']
    symbols = [cast(Optional[str], number_to_symbol.get(i)) for i in range(COMMON_ERROR_FIRST, COMMON_ERROR_LAST)]
    common_entries_list = [f'"{elm.lower()}"' if elm is not None else 'nullptr' for elm in symbols]
    common_entries = ''.join(f'    {elm},\n' for elm in common_entries_list)

    # DB specific entries
    def _gen_specific_entries(df_db: pd.DataFrame) -> str:
        return ''.join(SERVER_ERROR_TO_STRING_TEMPLATE_ENTRY.format(
            number=r.numbr,
            symbol_lower=r.symbol.lower()
        ) for r in df_db.itertuples())
    mysql_entries = _gen_specific_entries(df_mysql)
    mariadb_entries = _gen_specific_entries(df_mariadb)
     
    return SERVER_ERROR_TO_STRING_TEMPLATE.format(
        common_error_first=COMMON_ERROR_FIRST,
        common_entries=common_entries,
        mysql_entries=mysql_entries,
        mariadb_entries=mariadb_entries
    )


# Parse a header into a dataframe of (number, symbol) pairs
def parse_err_header(fname: Path) -> pd.DataFrame:
    with open(fname, 'rt') as f:
        lines = f.read().split('\n')
    v = [elm for elm in lines if elm.startswith('#define')]
    v = [elm.split(' ')[1:] for elm in v]
    df = pd.DataFrame(v, columns=['symbol', 'numbr'])
    df = df[~df['numbr'].isna()]
    df['numbr'] = df['numbr'].astype(int)
    df = df[df['numbr'] < SERVER_ERROR_LAST]
    # Discard pseudo error codes that some header have
    df = df[df['symbol'].map(lambda x: not (
        x.startswith('ER_ERROR_FIRST') or
        x.startswith('ER_ERROR_LAST') or
        x == 'ER_LAST_MYSQL_ERROR_MESSAGE' or
        x.startswith('ER_UNUSED') or
        x.endswith('__UNUSED')
    ))]
    return df


# MySQL 5.x and 8.x don't fully agree on error names. Some names have been
# removed, others have been added and others have been renamed. We merge
# both so the library can be used with both systems. In case of conflict, pick the 8.x name
# (they generally add a _UNUSED suffix for the codes they no longer use).
# Some symbols appear both in 5.x and 8.x but with different values - we pick the 8.x in
# case of conflict. 
def merge_mysql_errors(df_mysql5: pd.DataFrame, df_mysql8: pd.DataFrame) -> pd.DataFrame:
    def resolve_symbol(r):
        s5 = r['symbol_5']
        s8 = r['symbol_8']
        if not pd.isna(s5) and pd.isna(s8):
            symbol, dbver = s5, 5
        else:
            symbol, dbver = s8, 8
        return pd.Series(dict(numbr=r['numbr'], symbol=symbol, dbver=dbver))
    
    return df_mysql5 \
        .join(df_mysql8.set_index('numbr'), how='outer', on='numbr', lsuffix='_5', rsuffix='_8') \
        .apply(resolve_symbol, axis=1) \
        .sort_values(by='dbver') \
        .drop_duplicates(['symbol'], keep='last') \
        .drop(columns=['dbver'])


# Split between common and specific codes
def generate_error_ranges(df_mysql: pd.DataFrame, df_mariadb: pd.DataFrame) -> ErrorCodes:
    # Join
    joined = df_mysql.join(df_mariadb.set_index('numbr'), how='outer', on='numbr', lsuffix='_mysql', rsuffix='_mariadb')
    joined = joined[joined['numbr'] < COMMON_ERROR_LAST]

    # Common range
    res_common = joined[joined['symbol_mysql'] == joined['symbol_mariadb']]
    res_common = res_common.rename(columns={'symbol_mysql': 'symbol'}).drop(columns=['symbol_mariadb'])

    # Values in the common range that differ between mysql and mariadb
    joined_different = joined[joined['symbol_mysql'] != joined['symbol_mariadb']]
    res_mysql_1 = joined_different[joined_different['symbol_mysql'].notna()].rename(columns={'symbol_mysql': 'symbol'}).drop(columns=['symbol_mariadb'])
    res_mariadb_1 = joined_different[joined_different['symbol_mariadb'].notna()].rename(columns={'symbol_mariadb': 'symbol'}).drop(columns=['symbol_mysql'])

    # Values that are outside the common range
    res_mysql_2 = df_mysql[df_mysql['numbr'] >= COMMON_ERROR_LAST]
    res_mariadb_2 = df_mariadb[df_mariadb['numbr'] >= COMMON_ERROR_LAST]

    return ErrorCodes(
        common=res_common.sort_values(by='numbr'),
        mysql=pd.concat([res_mysql_1, res_mysql_2]).sort_values(by='numbr'),
        mariadb=pd.concat([res_mariadb_1, res_mariadb_2]).sort_values(by='numbr')
    )

# Does the full process. folder should contain the relevant headers
def parse_headers(folder: Path) -> ErrorCodes:
    df_mysql8_header = parse_err_header(folder.joinpath('mysql8.h'))
    df_mysql5_header = parse_err_header(folder.joinpath('mysql5.h'))
    df_mariadb_header = parse_err_header(folder.joinpath('mariadb.h'))
    df_mysql_header = merge_mysql_errors(df_mysql5_header, df_mysql8_header)
    return generate_error_ranges(df_mysql_header, df_mariadb_header)


# Writes a CSV file with the contents of headers, so we can keep track of what we did in the last Boost version
def write_csv(codes: ErrorCodes) -> None:
    df = pd.concat([
        codes.common.assign(category='common'),
        codes.mysql.assign(category='mysql'),
        codes.mariadb.assign(category='mariadb')
    ]).sort_values(by=['numbr', 'category'])
    df.to_csv(_CSV_PATH, index=False)


# Loads the CSV file from the previous version
def load_csv() -> ErrorCodes:
    df = pd.read_csv(_CSV_PATH)
    return ErrorCodes(
        common=df[df['category'] == 'common'].drop(columns=['category']),
        mysql=df[df['category'] == 'mysql'].drop(columns=['category']),
        mariadb=df[df['category'] == 'mariadb'].drop(columns=['category']),
    )


def _merge_new_codes_single(df_common: pd.DataFrame, df_old: pd.DataFrame, df_new: pd.DataFrame) -> pd.DataFrame:
    # Remove anything that's already present in the current headers
    temp = pd.concat([df_common, df_new, df_common]).drop_duplicates(subset='numbr', keep=False)
    temp = pd.concat([df_old, temp, df_old]).drop_duplicates(subset='numbr', keep=False)
    return pd.concat([df_old, temp])


def merge_new_codes(old: ErrorCodes, new: ErrorCodes) -> ErrorCodes:
    return ErrorCodes(
        common=old.common, # the common range never gets modified
        mysql=_merge_new_codes_single(old.common, old.mysql, new.mysql),
        mariadb=_merge_new_codes_single(old.common, old.mariadb, new.mariadb),
    )


# Actually perform the generation
def write_headers(codes: ErrorCodes) -> None:
    def header_path(p: List[str]) -> Path:
        return REPO_BASE.joinpath('include', 'boost', 'mysql', *p)

    # common_server_errc.hpp
    with open(header_path(['common_server_errc.hpp']), 'wt') as f:
        f.write(render_common_server_errc(codes.common))
    
    # mysql_server_errc.hpp
    with open(header_path(['mysql_server_errc.hpp']), 'wt') as f:
        f.write(render_server_specific_errc('mysql', codes.mysql))

    # mariadb_server_errc.hpp
    with open(header_path(['mariadb_server_errc.hpp']), 'wt') as f:
        f.write(render_server_specific_errc('mariadb', codes.mariadb))
    
    # detail/auxiliar/server_errc_strings.hpp
    with open(header_path(['impl', 'internal', 'error', 'server_error_to_string.ipp']), 'wt') as f:
        f.write(render_server_error_to_string(codes.common, codes.mysql, codes.mariadb))


# We need to run file_headers.py to set copyrights and headers
def invoke_file_headers() -> None:
    run(['python', str(REPO_BASE.joinpath('tools', 'scripts', 'file_headers.py'))])


def main():
    old_codes = load_csv()
    new_codes = parse_headers(REPO_BASE.joinpath('private', 'errors', '1.86'))
    codes = merge_new_codes(old_codes, new_codes)
    write_csv(codes)
    write_headers(codes)
    invoke_file_headers()

            
if __name__ == '__main__':
    main()
