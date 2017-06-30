#!/bin/sh

pdb=$(readlink -f $1)

if [ -r $pdb.gz ]; then
    gzip -d $pdb.gz
fi

cd $RWPLUS_HOME
./calRWplus $pdb | cut -d= -f2
