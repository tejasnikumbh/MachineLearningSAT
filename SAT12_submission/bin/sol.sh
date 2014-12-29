#!/bin/bash

# Top-level script for SOL based on RSat 2.01 

if [ "x$1" = "x" ]; then
  echo "c USAGE: sol.sh <input CNF>"
  exit 1
fi

path=`dirname $0` 
INPUT=$1;
shift

mono $path/sol.exe $INPUT 
