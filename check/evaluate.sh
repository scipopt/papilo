#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

export LANG=C
export LC_NUMERIC=C

FILE=$1

BASENAME=`basename $FILE .out`
DIR=`dirname $FILE`
OUTFILE=$DIR/$BASENAME.out
ERRFILE=$DIR/$BASENAME.err
RESFILE=$DIR/$BASENAME.res
TEXFILE=$DIR/$BASENAME.tex
PAVFILE=$DIR/$BASENAME.pav

# detect test set
TSTNAME=`echo $BASENAME | sed 's/check.\([a-zA-Z0-9_-]*\).*/\1/g'`

# detect test used solver
SOLVER=`echo $BASENAME | sed 's/check.\([a-zA-Z0-9_-]*\).\([a-zA-Z0-9_]*\).*/\2/g'`

echo "Testset " $TSTNAME
echo "Solver  " $SOLVER

if test -f testset/$TSTNAME.test
then
   TESTFILE=testset/$TSTNAME.test
else
   TESTFILE=""
fi

# call method to obtain solution file
# defines the following environment variable: SOLUFILE
. ./configuration_solufile.sh $TSTNAME

echo $OUTFILE

# we should not check for SOLVER = scip here, because check.awk needs also to be called for examples with other names
awk -f check.awk -v "TEXFILE=$TEXFILE" -v "PAVFILE=$PAVFILE" -v "ERRFILE=$ERRFILE" $AWKARGS $TESTFILE $SOLUFILE $OUTFILE | tee $RESFILE
