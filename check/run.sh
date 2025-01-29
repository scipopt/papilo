#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2020-2025 Zuse Institute Berlin (ZIB)                    *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# The script executing EXECNAME on one instance and producing the logfiles.
# Can be executed either locally or on a cluster node.
# Is to be invoked inside a 'check(_cluster)*.sh' script.

# absolut tolerance for checking linear constraints and objective value
LINTOL=1e-04
# absolut tolerance for checking integrality constraints
INTTOL=1e-04

TMPDIR="${CLIENTTMPDIR}/${USER}-tmpdir"

# check if tmp-path exists
if test ! -d "${TMPDIR}"
then
    mkdir "${TMPDIR}"
    echo "Creating directory ${TMPDIR} for temporary outfile"
fi

OUTFILE="${TMPDIR}${BASENAME}.out"
ERRFILE="${TMPDIR}${BASENAME}.err"
SOLFILE="${TMPDIR}${BASENAME}.sol"
DATFILE="${TMPDIR}${BASENAME}.dat"
TMPFILE="${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.tmp"

uname -a                            > "${OUTFILE}"
uname -a                            > "${ERRFILE}"

# function to copy back the results and delete temporary files
function cleanup {
    mv "${OUTFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.out"
    mv "${ERRFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.err"
    if [ -d "${ERRFILE}.rr" ]
    then
        mv -f "${ERRFILE}.rr" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.rr"
    fi
    # move a possible data file
    if [ -f "${DATFILE}" ] ;
    then
        mv "${DATFILE}" "${SOLVERPATH}/${OUTPUTDIR}/${BASENAME}.dat"
    fi
    rm -f "${TMPFILE}"
    rm -f "${SOLFILE}"
}

# ensure TMPFILE is deleted and results are copied when exiting (normally or due to abort/interrupt)
trap cleanup EXIT

# only wait for optimi to be mounted in run.sh if you are on an opt computer at zib
OPTHOST=$(uname -n | sed 's/.zib.de//g' | sed 's/portal//g' | tr -cd '[:alpha:]')

# check if the scripts runs a *.zib.de host
if $(hostname -f | grep -q zib.de) && $([[ "${OPTHOST}" == "opt" ]] || [[ "${OPTHOST}" == "optc" ]]);
then
    # access /optimi once to force a mount
    ls /nfs/optimi/QUOTAS >/dev/null 2>&1

    # check if /optimi is mounted
    MOUNTED=0

    # count number of fails and abort after 10 min to avoid an endless loop
    FAILED=0

    while [ "${MOUNTED}" -ne 1 ]
    do
        # stop if the system does not mount /optimi for ~10 minutes
        if [ "${FAILED}" -eq 600 ]
        then
            exit 1
        fi

        if [ -f /nfs/optimi/QUOTAS ] ;
        then
            MOUNTED=1
        else
            ((FAILED++))
            echo "/optimi is not mounted yet, waiting 1 second"
            sleep 1
        fi
    done
fi

echo                                      >> "${OUTFILE}"
if test $(uname) == Linux ; then   # -b does not work with top on macOS
    top -b -n 1 | head -n 15              >> "${OUTFILE}"
fi
echo                                      >> "${OUTFILE}"
echo "hard time limit: ${HARDTIMELIMIT}"  >> "${OUTFILE}"
echo "hard mem limit: ${HARDMEMLIMIT}"    >> "${OUTFILE}"
echo                                      >> "${OUTFILE}"
echo "SLURM jobID: ${SLURM_JOB_ID}"       >> "${OUTFILE}"
echo                                      >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="        >> "${OUTFILE}"
echo "@01 ${FILENAME} ==========="        >> "${ERRFILE}"
echo "-----------------------------"      >> "${OUTFILE}"
date                                      >> "${OUTFILE}"
date                                      >> "${ERRFILE}"
echo "-----------------------------"      >> "${OUTFILE}"
date +"@03 %s"                            >> "${OUTFILE}"
echo "@05 ${TIMELIMIT}"                   >> "${OUTFILE}"

#if we use a debugger command, we need to replace the errfile place holder by the actual err-file for logging
#and if we run on the cluster we want to use srun with CPU binding which is defined by the check_cluster script
EXECNAME="${EXECNAME/ERRFILE_PLACEHOLDER/${ERRFILE}}"
EXECNAME="${SRUN}${EXECNAME/RRTRACEFOLDER_PLACEHOLDER/${ERRFILE}}"

if [[ ${SOLVE_EXECUTABLE} == "off" ]]
then
    # FILENAME is exported in the calling script check_cluster.sh
    echo ">>> Executing: ${EXECNAME} ${PAPILO_OPT_COMMAND} -f ${FILENAME} -p ${SETFILEPAPILO} -s ${SETFILESCIP} --tlim ${TIMELIMIT} --presolve.randomseed=${SEED}"
    eval "${EXECNAME} ${PAPILO_OPT_COMMAND} -p ${SETFILEPAPILO} -s ${SETFILESCIP} -f ${FILENAME} --tlim ${TIMELIMIT} --presolve.randomseed=${SEED}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
else
    if [[ ${SKIP_PRESOLVE} == "false" ]]
    then
        PRESOLVED_FILENAME="${FILENAME[@]//.opb/_presolved.opb}"
        echo ">>> Executing: ${EXECNAME} presolve -f ${FILENAME} -p ${SETFILEPAPILO} -s ${SETFILESCIP} --tlim ${TIMELIMIT} --presolve.randomseed=${SEED} -r ${PRESOLVED_FILENAME}"
        eval "${EXECNAME} presolve -p ${SETFILEPAPILO} -s ${SETFILESCIP} -f ${FILENAME} --tlim ${TIMELIMIT} --presolve.randomseed=${SEED} -r ${PRESOLVED_FILENAME}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
    else
        echo "skipped presolving"
        PRESOLVED_FILENAME="${FILENAME}"
    fi

    if [[ ${SOLVE_EXECUTABLE} =~ $"scip" ]]
    then
        echo ">>> Executing SCIP ${SOLVE_EXECUTABLE} -f ${PRESOLVED_FILENAME}"
        eval "${SOLVE_EXECUTABLE} -f ${PRESOLVED_FILENAME}" -s ${SETFILESCIP} 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
    elif [[ ${SOLVE_EXECUTABLE} =~ $"roundingsat" ]]
    then
        echo ">>> Executing ROUNDINGSAT ${SOLVE_EXECUTABLE} ${PRESOLVED_FILENAME}"
        eval "${SOLVE_EXECUTABLE} ${PRESOLVED_FILENAME}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
    elif [[ ${SOLVE_EXECUTABLE} =~ $"soplex" ]]
    then
        echo ">>> Executing SOPLEX ${SOLVE_EXECUTABLE} ${PRESOLVED_FILENAME}"
        eval "${SOLVE_EXECUTABLE} ${PRESOLVED_FILENAME}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
    elif [[ ${SOLVE_EXECUTABLE} =~ .*"sat4j-pb.".* ]]
    then
        echo ">>> Executing Sat4j java -jar ${SOLVE_EXECUTABLE} ${PAPILO_OPT_COMMAND} ${TIMELIMIT} ${PRESOLVED_FILENAME}"
        eval "java -jar ${SOLVE_EXECUTABLE} ${PAPILO_OPT_COMMAND} ${TIMELIMIT} ${PRESOLVED_FILENAME}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
    else
      echo "Unknown solver: ${SOLVE_EXECUTABLE}"
    fi
    if [[ ${SKIP_PRESOLVE} == "false" ]]
    then
      rm ${PRESOLVED_FILENAME}
    fi
fi

# TODO extract these variables to be configurable from make file
VERIPB_PATH="/home/alexander/git_repositories/veripb-dev/venv/bin/activate"
VERIPB=true
if ${VERIPB} && ( test -z "${VERIPB_PATH}" || test -f ${VERIPB_PATH} )
then
  if [[ ${FILENAME} == *.opb.gz ]]
  then
      eval gzip -dkf ${FILENAME}
      PBP_FILENAME=${FILENAME%.opb.gz}.pbp
      if test -f ${PBP_FILENAME}
      then
        UNZIPPED_FILENAME=${FILENAME%.opb.gz}.opb
        echo ">>> Executing: veripb --stats --forceCheckDeletion ${UNZIPPED_FILENAME} ${PBP_FILENAME}"
        if ! test -z "${VERIPB_PATH}"
        then
          eval "source ${VERIPB_PATH}"
        fi
        eval "veripb --stats --forceCheckDeletion ${UNZIPPED_FILENAME} ${PBP_FILENAME}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
        eval "rm -r ${UNZIPPED_FILENAME}"
      fi
  elif [[ ${FILENAME} == *.opb ]]
  then
      PBP_FILENAME=${FILENAME%.opb}.pbp
      if test -f ${PBP_FILENAME}
      then
        echo ">>> Executing: veripb --stats --forceCheckDeletion ${FILENAME} ${PBP_FILENAME}"
        if ! test -z "${VERIPB_PATH}"
        then
          eval "source ${VERIPB_PATH}"
        fi
        eval "veripb --stats --forceCheckDeletion ${FILENAME} ${PBP_FILENAME}" 2>> "${ERRFILE}" | tee -a "${OUTFILE}"
      fi
    else
      echo "$FILENAME has to end with .opb or .opb.gz to run Veripb"
  fi
fi
echo "----------------------------------------------------------"

retcode=${PIPESTATUS[0]}
if test "${retcode}" != 0
then
    echo "${EXECNAME} returned with error code ${retcode}." >> "${ERRFILE}"
fi

if test -e "${SOLFILE}"
then
    # translate SCIP solution format into format for solution checker. The
    # SOLFILE format is a very simple format where in each line we have a
    # <variable, value> pair, separated by spaces.  A variable name of
    # =obj= is used to store the objective value of the solution, as
    # computed by the solver. A variable name of =infeas= can be used to
    # indicate that an instance is infeasible.
    sed ' /solution status:/d;
    s/objective value:/=obj=/g;
    s/infinity/1e+20/g;
    s/no solution available//g' "${SOLFILE}" > "${TMPFILE}"
    mv "${TMPFILE}" "${SOLFILE}"

fi

date +"@04 %s"                        >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${OUTFILE}"
echo "-----------------------------"  >> "${OUTFILE}"
date                                  >> "${ERRFILE}"
echo                                  >> "${OUTFILE}"
echo "=ready="                        >> "${OUTFILE}"
