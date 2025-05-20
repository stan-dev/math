#!/usr/bin/python3
#
# Copyright (c) 2019-2024 Ruben Perez Hidalgo (rubenperez038 at gmail dot com)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#

import argparse
from subprocess import PIPE, STDOUT, Popen
from contextlib import contextmanager
import re
import os
import socket
import struct

_is_win = os.name == 'nt'


# Returns the port the server is listening at
def _parse_server_start_line(line: str) -> int:
    m = re.match(r'Server listening at 0\.0\.0\.0:([0-9]+)', line)
    if m is None:
        raise RuntimeError('Unexpected server start line')
    return int(m.group(1))


@contextmanager
def _launch_server(exe: str, host: str):
    # Launch server and let it choose a free port for us.
    # This prevents port clashes during b2 parallel test runs
    server = Popen([exe, 'example_user', 'example_password', host, '0'], stdout=PIPE, stderr=STDOUT)
    assert server.stdout is not None
    with server:
        try:
            # Wait until the server is ready
            ready_line = server.stdout.readline().decode()
            print(ready_line, end='', flush=True)
            if ready_line.startswith('Sorry'): # C++ standard unsupported, skip the test
                exit(0)
            yield _parse_server_start_line(ready_line)
        finally:
            print('Terminating server...', flush=True)

            # In Windows, there is no sane way to cleanly terminate the process.
            # Sending a Ctrl-C terminates all process attached to the console (including ourselves
            # and any parent test runner). Running the process in a separate terminal doesn't allow
            # access to stdout, which is problematic, too.
            if _is_win:
                # kill is an alias for TerminateProcess with the given exit code
                os.kill(server.pid, 9999)
            else:
                # Send SIGTERM
                server.terminate()

            # Print any output the process generated
            print('Server stdout: \n', server.stdout.read().decode(), flush=True)
    
    # Verify that it exited gracefully
    if (_is_win and server.returncode != 9999) or (not _is_win and server.returncode):
        raise RuntimeError('Server did not exit cleanly. retcode={}'.format(server.returncode))


class _Runner:
    def __init__(self, port: int) -> None:
        self._port = port
    
    def _connect(self) -> socket.socket:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.connect(('127.0.0.1', self._port))
        return sock


    def _query_employee(self, employee_id: int) -> str:
        # Open a connection
        sock = self._connect()

        # Send the request
        sock.send(struct.pack('>Q', employee_id))
        
        # Receive the response. It should always fit in a single TCP segment
        # for the values we have in CI
        res = sock.recv(4096).decode()
        assert len(res) > 0
        return res


    def _generate_error(self) -> None:
        # Open a connection
        sock = self._connect()

        # Send an incomplete message
        sock.send(b'abc')
        sock.close()
    

    def run(self, test_errors: bool) -> None:
        # Generate an error first. The server should not terminate
        if test_errors:
            self._generate_error()
        assert self._query_employee(1) != 'NOT_FOUND'
        value = self._query_employee(0xffffffff)
        assert value == 'NOT_FOUND', 'Value is: {}'.format(value)


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('host')
    parser.add_argument('--test-errors', action='store_true')
    args = parser.parse_args()

    # Launch the server
    with _launch_server(args.executable, args.host) as listening_port:
    # Run the tests
        _Runner(listening_port).run(args.test_errors)


if __name__ == '__main__':
    main()
