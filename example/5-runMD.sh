#!/bin/bash

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
RUN_SCRIPT_PATH=$SHELL_FOLDER/../runMolecularDynamics/run_folder/Gm12878/chr1/run00
RUN_SCRIPT=$RUN_SCRIPT_PATH/run.sh

echo '''
****** Run Simulation ******'''
echo '
>>>> This script is to run the LAMMPS simulation for chromosome 1, GM12878.'

cd $RUN_SCRIPT_PATH
./run.sh
