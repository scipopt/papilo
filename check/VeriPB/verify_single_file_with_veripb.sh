#!/bin/bash -ex
echo "${FILENAME}"
instance=${FILENAME%.pbp}.opb
#  sbatch -c "veripb --trace --useColor ${instance} ${f} > ${out}"
echo ">>> Executing: veripb --trace --useColor ${instance} ${FILENAME}"
source "/scratch/opt/ahoen/VeriPB/z1venv/bin/activate"
eval "veripb --trace --useColor ${instance} ${FILENAME}"
