#!/bin/bash -ex
for f in /data/optimi/optimi/kombadon/IP/MIPLIB01/dec/miplib2017/*.pbp;
do
  export FILENAME="${f}"
  echo "${FILENAME}"
  out=${FILENAME%.pbp}.out
  sbatch --constraint="Gold5122" --partition="opt_int" -A "optimi_integer" --output=${out} --cpu-freq=highm1  --export=FILENAME=$FILENAME verify_single_file_with_veripb.sh
done
