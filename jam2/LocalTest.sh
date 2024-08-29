#!/bin/bash
#./jam -f input/TestCascadeDecay.inp
./jam -f $1
root4star -b -l -q makeTree.C+\(\"phase.txt\"\)
