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
        g++ -std=c++11 -pthread $1 -I./ -I$SROOT/include -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include -I$GOLEMSPACE/include -I$GOLEMSPACE/local/lib64 -L$SROOT/lib -L$SROOT/lib64 -L/usr/include/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/lib -L$GOLEMSPACE/local/lib  -lSQuIDS -lnuSQuIDS -lhdf5_hl -lhdf5 -lgsl -lgslcblas -o $2
    fi
fi


