#!/bin/bash

DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
ROOTDIR=$DIR/../..
DOCSDIR=$DIR/../../edog-docs

if [ -d $DOCSDIR/.git ]; then
  echo "All good, found edog-docs with git repo."
  echo "Pulling changes in edog-docs (note: not pulling in the current folder)."
  cd $DOCSDIR
  if ! git reset --hard; then
      echo "Could not reset"
      exit 1
  fi
  if ! git pull; then
      echo "Could not pull"
      exit 1
  fi
elif [ -d $DOCSDIR ]; then
  echo "Directory $DOCSDIR exists, but is not a git repo. Please delete the directory first."
  exit 1
else
  cd $ROOTDIR
  if ! git clone -b gh-pages git@github.com:miladh/extendedDOG.git edog-docs; then
      echo "Could not clone"
      exit 1
  fi
fi

cd $DIR
echo Running $DOCSDIR/qdoc edog.qdocconf
LD_LIBRARY_PATH=$DOCSDIR $DOCSDIR/qdoc edog.qdocconf
