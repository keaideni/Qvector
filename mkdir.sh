#!/bin/sh

mkdir block
mkdir result
mkdir trunc
cat>Parameter<<EOF
nmax= 6 
D= 100
LatticeSize= 100
omega0= 1
omegaq= 1
gr= 0.1
gcr= 0.1
Jr= 0.1
Jcr= 0.1
EOF
