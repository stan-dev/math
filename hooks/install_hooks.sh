#!/bin/bash

for i in `find $(dirname "$0") -type f | grep -v 'install_hooks.sh'`; do
    dest=`git rev-parse --git-dir`/hooks
    echo "linking $i to $dest"
    chmod +x $i
    ln -s `pwd`/$i $dest/$(basename $i)
done
