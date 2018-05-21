#!/bin/csh -fx

#### Make makdep executable
gmake

#### Copy it back on directory to main SRC directory so it can act on SRC subdirectories
cp makdep ..

