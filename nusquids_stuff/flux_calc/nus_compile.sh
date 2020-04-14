#!/bin/bash

# Ben Smithers

# short script to make it easier to compile individual nusquids scripts 
# arg1 - file to compile
# arg2 - output binary file

if [ -z "$1" ]
then
    echo "Nothing given to compile!"
    exit 1
else
    if [ -z "$2" ]
    then 
        echo "Specify output file!"
        exit 1
    else
        echo "Building $1 into $2"
        g++ -std=c++11 $1 -I./ -I$SROOT/include -I$GOLEMSPACE/include -L$SROOT/lib -L$SROOT/lib64 -L$GOLEMSPACE/local/lib -I$GOLEMSPACE/local/lib64 -lSQuIDS -lnuSQuIDS -lhdf5_hl -lhdf5 -lgsl -lgslcblas -o $2
    fi
fi


