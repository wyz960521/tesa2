#/bin/bash

# Ensure that the current directory is the TESA root directory

# TESA file compilation
cd src
make clean && make
cd ..

# Running TESA using an input file without sequencing coverage
cd src
./tesa -i ../example/test.fasta -l 14
./tesa -i ../example/test.fasta -L 14 -U 16
