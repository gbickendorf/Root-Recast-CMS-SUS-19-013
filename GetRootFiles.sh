#!/bin/bash
if [ ${#*} -ne 1 ]; then
    echo "usage: `basename $0` basedir"
    exit 1
fi
find $(realpath $1) -iname "*.root" >ROOTFILES.txt
wc -l ROOTFILES.txt
