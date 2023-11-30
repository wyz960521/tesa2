#/bin/bash

# Ensure that the current directory is the tesa root directory

chmod +x preprocess.sh
cd example
../preprocess.sh test reference.fa reference.fa.fai test_out
../src/tesa test_out.tesa
