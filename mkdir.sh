#!/bin/sh

mkdir block
mkdir result
mkdir trunc
cat>Parameter<<EOF
nmax= 6 
D= 200
LatticeSize= 10
omega0= 1
omegaq= 1
gr= 2
gcr= 0
Jr= 0
Jcr= 0
EOF
