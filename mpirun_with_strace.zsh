#!/usr/bin/zsh

#!/bin/bash
mpirun -np 8 strace_proc.zsh pol/1PGB.pol pot/1PGB.pot provatot.par
