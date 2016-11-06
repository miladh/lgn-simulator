#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
ROOTDIR=$DIR/../..
DOCSDIR=$DIR/../../lgn-simulator-docs



cd $DIR
echo Running $DOCSDIR/qdoc lgn-simulator.qdocconf
LD_LIBRARY_PATH=$DOCSDIR $DOCSDIR/qdoc lgn-simulator.qdocconf
