#!/bin/bash

BASE_DIR=/storage/gpfs_data/limadou/puccetti/OLD_WORK/ML_Project/ml_project/cluster
EXE=${BASE_DIR}/src/stabilitytime
OUTPUT_DIR=${BASE_DIR}/output/data

JOB_ID=$1

# Grid of tunes: 15 values in [0.10, 0.45] step 0.025
NU_VALUES=(0.100 0.125 0.150 0.175 0.200 0.225 0.250 0.275 0.300 0.325 0.350 0.375 0.400 0.425 0.450)

# N values (log-spaced)
N_VALUES=(100000 300000 1000000 3000000 10000000 30000000 100000000)

N_NU=${#NU_VALUES[@]}   # 15
N_N=${#N_VALUES[@]}     # 7
# Total jobs: 15 * 15 * 7 = 1575

# Decode JOB_ID -> (i_nux, i_nuy, i_N)
TUNE_PAIR=$(( JOB_ID / N_N ))         # 0..224
I_N=$(( JOB_ID % N_N ))               # 0..6
I_NUX=$(( TUNE_PAIR / N_NU ))         # 0..14
I_NUY=$(( TUNE_PAIR % N_NU ))         # 0..14

NU_X=${NU_VALUES[$I_NUX]}
NU_Y=${NU_VALUES[$I_NUY]}
N=${N_VALUES[$I_N]}

EPS=0
N_STEPS=300
LATTICE_LENGTH=0.45

echo "Job ${JOB_ID}: nu_x=${NU_X} nu_y=${NU_Y} eps=${EPS} N=${N}"

cd ${OUTPUT_DIR}
${EXE} ${NU_X} ${NU_Y} ${EPS} ${N_STEPS} ${LATTICE_LENGTH} ${N}
