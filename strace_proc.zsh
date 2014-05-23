#!/bin/zsh
/usr/bin/strace -T -ttt -o /tmp/strace.out.$$ ./montegrappa $@
