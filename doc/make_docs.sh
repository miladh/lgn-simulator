#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
ROOTDIR=$DIR/../..
DOCSDIR=$DIR/../../edog-docs



cd $DIR
echo Running $DOCSDIR/qdoc edog.qdocconf
LD_LIBRARY_PATH=$DOCSDIR $DOCSDIR/qdoc edog.qdocconf
