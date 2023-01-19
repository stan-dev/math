import subprocess
from pathlib import Path

def get_dependencies(file: Path) -> set[str]:
    file_dot_d = file.with_suffix('.d')
    subprocess.run(['make', str(file_dot_d)], stdout=subprocess.DEVNULL)
    contents = file_dot_d.read_text()
    file_dot_d.unlink()

    ## TODO include self as a 'dependency'
    return set(filter(lambda s: 'stan/math' in s and s.endswith(".hpp"), map(str.strip, contents.split(' '))))

def get_changed() -> set[str]:
    changed_files = subprocess.run(
        ["git", "diff", "--name-only", "origin/develop...HEAD"], text=True, capture_output=True
    ).stdout.splitlines()
    return set(changed_files)

# recursively ls test/prob/*
# test each file with get_dependencies
# if file has changed deps, run runTests on $file.replace("._test.hpp", "_*_test.hpp")

changed = get_changed()
tests_to_run = []

distribution_tests = Path("test", "prob")

for dist in distribution_tests.iterdir():
    if not dist.is_dir():
        continue
    for test in dist.iterdir():
        deps = get_dependencies(test)
        if changed & deps:
            tests_to_run.append(str(test).replace("_test.hpp", "_0*_test.cpp"))

print('\n'.join(tests_to_run))
