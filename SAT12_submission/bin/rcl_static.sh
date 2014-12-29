#!/bin/bash

# Top-level script for RSat 2.01 

if [ "x$1" = "x" ]; then
  echo "c USAGE: SatELiteGTI <input CNF>"
  exit 1
fi

path=`dirname $0` 
TMP=/tmp/$$_temp_many  #set this to the location of temporary files
SE=$path/SatELite_release_RCL           #set this to the executable of SatELite
RS=$path/"cdcl"               #set this to the executable of RSat
INPUT=$1;
shift

$SE $INPUT $TMP.cnf $TMP.vmap $TMP.elim
X=$?


if [ $X == 20 ]
then
    echo "s UNSATISFIABLE"
    rm -f $TMP.cnf $TMP.vmap $TMP.elim $TMP.result
    exit 20
fi

if [ $X == 0 ]; then
  #SatElite terminated correctly
    $RS < $TMP.cnf  > $TMP.result 2>/dev/null
    X=$?
    
    if [ $X == 20 ]; then
    #RSat must not print out result!
	echo "s UNSATISFIABLE"
	rm -f $TMP.cnf $TMP.vmap $TMP.elim $TMP.result
	exit 20
    #Don't call SatElite for model extension.
    elif [ $X != 10 ]; then
    #timeout/unknown, nothing to do, just clean up and exit.
	rm -f $TMP.cnf $TMP.vmap $TMP.elim $TMP.result
	exit $X
    fi 
  #SATISFIABLE, call SatElite for model extension
    echo SAT > $TMP.result1
    cat $TMP.result | grep v | tail -c+2 >> $TMP.result1
    mv $TMP.result1 $TMP.result
    cat $TMP.result
    $SE +ext $INPUT $TMP.result $TMP.vmap $TMP.elim
    X=$?
    rm -f $TMP.cnf $TMP.vmap $TMP.elim $TMP.result $TMP.model*
else        
    $RS < $INPUT  #but we must force rsat to print out result here!!!
    X=$?
fi

rm -f $TMP.cnf $TMP.vmap $TMP.elim $TMP.result
exit $X
