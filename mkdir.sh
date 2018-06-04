#!/bin/sh

mkdir block
mkdir result
mkdir trunc
cat>Parameter<<EOF
nmax= 6
D= 20
LatticeSize= 10
gr= 1
gcr= 1
Jr= 100
Jcr= 0
EOF
