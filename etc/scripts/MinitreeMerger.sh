#!/bin/bash

usage() {
    echo "MinitreeMerger.sh dir 0/1 [default=0]"
    echo "  - 2d arg = 1: remove *txt and job_*root files from the dir"
}

if [ ! $# -ge 1 ]; then
    usage
fi

dir=$1
remove=0
if [ $# -ge 2 ]; then
 remove=$2
fi


dirtmp=`pwd`
echo "Merging data:"
cd $dir/data
rm -f MiniTreeData.root
hadd MiniTreeData.root *.root
if [ $remove -eq 1 ]; then
    rm -f job_*root
    rm -f job_*root*txt
fi

cd $dirtmp