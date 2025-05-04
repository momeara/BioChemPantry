#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N search-sea

### This is a SGE submission script to screen targets against a library.
###
###   Directories and files used by thi script
###     <run_base>/
###        inputs/
###          <target_id>.csv   # with column headers (compound, smiles)
###        logs/
###        outputs/
###        sea-wrapper.sh
###     /scratch/$(whoami)/
###        <job_id>/<job_id>/<task_name>
###
###     <library_fname> # .sea library
###  
###   To submit job
###     cd <run_base>/logs
###     qsub -t 1-<n_targets> ../sea-wrapper.sh <run_base> <library_fname>
###     



set -e

PERSIST=$1

LIB_FULL=$2
LIB=$(basename $LIB_FULL)


INPUT_FILES=$PERSIST/inputs
COMPLETE=$PERSIST/outputs
TASK_INPUT=$( ls $INPUT_FILES | sed -n ${SGE_TASK_ID}p )
TASK_NAME=$( basename ${TASK_INPUT%.*} )
TASK_FILE=$INPUT_FILES/$TASK_INPUT

SCRATCH_DIR=/scratch
if [ ! -d $SCRATCH_DIR ]; then
    SCRATCH_DIR=/tmp
fi

TASK_DIR=$SCRATCH_DIR/$( whoami )/$JOB_ID/$TASK_NAME

echo "SGE_TASK_ID:" $SGE_TASK_ID 1>&2
echo "Source file:" $TASK_FILE 1>&2
echo "Run dir:" $( hostname ):$TASK_DIR 1>&2
echo "SEA library:" $LIB_FULL 1>&2

echo "source ~/init_ml.sh"
source /mnt/nfs/home/momeara/init_ml.sh 1>&2

mkdir -pv $TASK_DIR 1>&2
pushd $TASK_DIR 1>&2

cp $LIB_FULL . 
SeaShell.py \
    batch \
    --library $LIB \
    --disable-precache \
    --skip-standardization \
    --target-name $TASK_NAME \
    $TASK_FILE \
    $TASK_INPUT.out.csv 1>&2
/bin/rm -f $LIB*

popd
mkdir -pv $COMPLETE 1>&2
mv -v $TASK_DIR $COMPLETE/$TASK_NAME 1>&2
rm -rvf $TASK_DIR 1>&2
