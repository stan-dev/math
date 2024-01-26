#!/usr/bin/python3
#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import requests
from bs4 import BeautifulSoup
import os
from os import path

REPO_BASE = path.abspath(path.join(path.dirname(__file__), '..', '..'))
DOC_PATH = path.join(REPO_BASE, 'doc', 'html')

def list_doc_files():
    all_files = []
    for base_dir, _, files in os.walk(DOC_PATH):
        all_files += [path.join(base_dir, f) for f in files if f.endswith('.html')]
    return all_files

def get_href(elm, current_file):
    try:
        res = elm['href']
    except KeyError:
        return None
    if res.startswith('http://') or res.startswith('https://'):
        if '#error_er_' in res:
            return res.split('#error_er_')[0]
        else:
            return res
    else:
        curdir = path.dirname(current_file)
        return path.realpath(path.join(curdir, res.split('#')[0]))

def extract_links():
    external_links = {}
    internal_links = {}
    
    for fname in list_doc_files():
        with open(fname, 'rt') as f:
            html_doc = f.read()
        soup = BeautifulSoup(html_doc, 'html.parser')
        links = [get_href(elm, fname) for elm in soup.find_all('a')]
        internal_links.update({ elm: fname for elm in links if elm is not None and elm.startswith('/')})
        external_links.update({ elm: fname for elm in links if elm is not None and \
                              (elm.startswith('http://') or elm.startswith('https://'))})
        
    return (external_links, internal_links)

def check_external_links(links):
    s = requests.Session()
    for url in sorted(links.keys()):
        print('Checking ', url)
        response = s.head(url, allow_redirects=True)
        if response.status_code != 200:
            print('  ++++ {} response code: {}'.format(url, response.status_code))
            
def check_internal_links(links):
    for target, link_file in links.items():
        if not path.exists(target):
            print('  ++++ Link {} in file {} does not exist'.format(target, link_file))
            
def main():
    external, internal = extract_links()
    check_external_links(external)
    check_internal_links(internal)
    
if __name__ == '__main__':
    main()
