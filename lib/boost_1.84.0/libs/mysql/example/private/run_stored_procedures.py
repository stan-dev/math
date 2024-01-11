#!/usr/bin/python3
#
# Copyright (c) 2019-2023 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

from subprocess import run, PIPE
import argparse
import re


def _parse_order_id(output: str) -> str:
    res = re.search(r'Order: id=([0-9]*)', output)
    assert res is not None
    return res.group(1)


def _parse_line_item_id(output: str) -> str:
    res = re.search(r'Created line item: id=([0-9]*)', output)
    assert res is not None
    return res.group(1)


class Runner:
    def __init__(self, exe: str, host: str) -> None:
        self._exe = exe
        self._host = host
    
    def run(self, subcmd: str, *args: str) -> str:
        cmdline = [self._exe, 'orders_user', 'orders_password', self._host, subcmd, *args]
        print(' + ', cmdline)
        res = run(cmdline, check=True, stdout=PIPE)
        print(res.stdout.decode())
        return res.stdout.decode()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('host')
    args = parser.parse_args()

    runner = Runner(args.executable, args.host)

    # Some examples require C++14 and produce an apology for C++11 compilers
    output = runner.run('get-products', 'feast')
    if "your compiler doesn't have the required capabilities to run this example" in output:
        print('Example reported unsupported compiler')
        exit(0)
    
    order_id = _parse_order_id(runner.run('create-order'))
    runner.run('get-orders')
    line_item_id = _parse_line_item_id(runner.run('add-line-item', order_id, '1', '5'))
    runner.run('add-line-item', order_id, '2', '2')
    runner.run('add-line-item', order_id, '3', '1')
    runner.run('remove-line-item', line_item_id)
    runner.run('get-order', order_id)
    runner.run('checkout-order', order_id)
    runner.run('complete-order', order_id)
    runner.run('get-orders')


if __name__ == '__main__':
    main()
