#!/bin/bash
if [ ${#*} -ne 1 ]; then
    echo "usage: `basename $0` basedir"
    exit 1
fi
find $1 -iname "*.root" >ROOTFILES.txt
