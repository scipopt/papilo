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

# if RBCLI_TAG is set, then it will be passed as tags to rbcli
export LANG=C
export LC_NUMERIC=C

REMOVE=0
UPLOAD=0
APPEND=0
EXPIRE=0
AWKARGS=""
FILES=""

for i in $@
do
    if test ! -e "${i}"
    then
        if test "${i}" = "-r"
        then
            REMOVE=1
        elif test "${i}" = "-U"
        then
            UPLOAD=1
        elif test "${i}" = "-E"
        then
            UPLOAD=1
            EXPIRE=1
        elif test "${i}" = "-R"
        then
            REMOVE=1
            UPLOAD=1
        elif test "${i}" = "-T"
        then
            REMOVE=1
            UPLOAD=1
            EXPIRE=1
        else
            AWKARGS="${AWKARGS} ${i}"
        fi
    else
        FILES="${FILES} ${i}"
    fi
done

for FILE in ${FILES}
do

    DIR=$(dirname ${FILE})
    EVALFILE=$(basename ${FILE} .eval)
    EVALFILE=$(basename ${EVALFILE} .out)

    OUTFILE="${DIR}/${EVALFILE}.out"
    ERRFILE="${DIR}/${EVALFILE}.err"
    SETFILESCIP="${DIR}/${EVALFILE}.scip.set"
    SETFILEPAPILO="${DIR}/${EVALFILE}.papilo.set"
    METAFILE="${DIR}/${EVALFILE}.meta"

    # check if the eval file exists; if this is the case construct the overall solution files
    if test -e "${DIR}/${EVALFILE}.eval"
    then
        # in case an output file exists, copy it away to save the results
        DATEINT=$(date +"%s")
        if test -e "${OUTFILE}"
        then
           cp "${OUTFILE}" "${OUTFILE}.old-${DATEINT}"
        fi
        if test -e "${ERRFILE}"
        then
           cp "${ERRFILE}" "${ERRFILE}.old-${DATEINT}"
        fi

        echo > "${OUTFILE}"
        echo > "${ERRFILE}"
        echo ""
        echo "create overall output and error file for ${EVALFILE}"

        # check first if all out and err files exist for this eval-file.
        NMISSING=0
        for i in $(cat "${DIR}/${EVALFILE}.eval") DONE
        do
            if test "${i}" = "DONE"
            then
                break
            fi

            for extension in out err
            do
                FILE="${i}.${extension}"
                if ! test -e "${FILE}"
                then
                    echo Missing ${FILE}
                    ((NMISSING++))
                fi
            done
        done

        if [ "${NMISSING}" -gt 0 -a "${REMOVE}" -eq 1 ]
        then
            echo "Exiting because ${NMISSING} out/err file $([ ${NMISSING} -gt 1 ] && echo "s are" || echo " is" ) missing, please rerun without the REMOVE flag"
            exit
        fi

        for i in $(cat "${DIR}/${EVALFILE}.eval") DONE
        do
            if test "${i}" = "DONE"
            then
                break
            fi

            FILE="${i}.out"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${OUTFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            else
                echo @01 "${FILE}" ==MISSING==  >> "${OUTFILE}"
                echo                        >> "${OUTFILE}"
            fi

            FILE="${i}.err"
            if test -e "${FILE}"
            then
                cat "${FILE}" >> "${ERRFILE}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            else
                echo @01 "${FILE}" ==MISSING==  >> "${ERRFILE}"
                echo                        >> "${ERRFILE}"
            fi

            FILE="${i}.scip.set"
            if test -e "${FILE}"
            then
                cp "${FILE}" "${SETFILESCIP}"
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi

            FILE="${i}.tmp"
            if test -e "${FILE}"
            then
                if test "${REMOVE}" = "1"
                then
                    rm -f "${FILE}"
                fi
            fi
        done

        if test "${REMOVE}" = "1"
        then
            rm -f "${DIR}/${EVALFILE}.eval"
        fi
    fi

    # check if the out file exists
    if test -e "${OUTFILE}"
    then
        # upload results to rubberband.zib.de
        if test "${UPLOAD}" = "1"
        then
            if test "${EXPIRE}" = "1"
            then
                RB_EXP_DATE=$(date '+%Y-%b-%d' -d "+6 weeks")
                rbcli -e "${RB_EXP_DATE}" up "${OUTFILE}" "${ERRFILE}" "${SETFILESCIP}" "${METAFILE}"
            else
                if test -z "${RBCLI_TAG}"
                then
                    rbcli up "${OUTFILE}" "${ERRFILE}" "${SETFILESCIP}" "${METAFILE}"
                else
                    rbcli --tags "${RBCLI_TAG}" up "${OUTFILE}" "${ERRFILE}" "${SETFILESCIP}" "${METAFILE}"
                fi
            fi
        fi
    fi
done
