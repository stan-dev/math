"""
Functions and constants shared by runTests.py, runChecks.py and runClangd.py.
"""

from __future__ import annotations

import glob
import os
import shlex
import subprocess
import sys
import time
from typing import NoReturn

# The repository root, which holds this file.
ROOT = os.path.dirname(os.path.realpath(__file__))

winsfx = ".exe"
testsfx = "_test.cpp"


def stopErr(msg: str, returncode: int) -> NoReturn:
    """Report an error message to stderr and exit with a given code."""
    sys.stderr.write(f"{msg}\n")
    sys.stderr.write(f"exit now ({time.strftime('%x %X %Z')})\n")
    sys.exit(returncode)


def run_command(
    command: str | list[str],
    *,
    capture: bool = False,
    check: bool = True,
    echo: bool = False,
) -> subprocess.CompletedProcess[str]:
    """
    Run command and return the finished process. A string runs through the
    shell; a list runs the program directly.
    @param capture: keep stdout and stderr as text on the result instead of
        passing them through to the terminal
    @param check: stop the script with stopErr if the command fails or the
        program cannot be started. Without it, a program that cannot be
        started returns returncode 127, as a shell would.
    @param echo: print a separator line and the command before running it
    """
    shown = command if isinstance(command, str) else shlex.join(command)
    if echo:
        print("------------------------------------------------------------")
        print(shown)
    try:
        proc = subprocess.run(
            command,
            shell=isinstance(command, str),
            capture_output=capture,
            text=True,
            check=False,
        )
    except OSError as e:
        if check:
            stopErr(f"{shown} failed: {e}", 127)
        return subprocess.CompletedProcess(command, 127, "", str(e))
    if check and proc.returncode != 0:
        details = f"\n{proc.stderr}" if capture and proc.stderr else ""
        stopErr(f"{shown} failed{details}", proc.returncode)
    return proc


def files_in_folder(folder: str, suffix: str = "") -> list[str]:
    """Returns a list of files ending in suffix in the folder and all
    its subfolders recursively. The folder can be written with
    wildcards as with the Unix find command.
    """
    files: list[str] = []
    for f in glob.glob(folder):
        if os.path.isdir(f):
            files.extend(files_in_folder(f + os.sep + "**", suffix))
        elif f.endswith(suffix):
            files.append(f)
    return files


def test_files_in_folder(folder: str) -> list[str]:
    """Returns a list of test files (*_test.cpp) in the folder and all
    its subfolders recursively. The folder can be written with
    wildcards as with the Unix find command.
    """
    return files_in_folder(folder, testsfx)
