#!/bin/bash

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
MAIN_SCRIPT_PATH=$SHELL_FOLDER/../runMolecularDynamics/inputFiles/epig_input/
MAIN_SCRIPT_CS_PATH=$MAIN_SCRIPT_PATH/chromStates/
MAIN_SCRIPT_CTCF_PATH=$MAIN_SCRIPT_PATH/ctcfSites/

echo '''
****** Extract Epigenomics Input ******'''
echo '
>>>> This script is to extract the model input of chromatin states and CTCF-binding sites for chromosome 1, GM12878.'

cd $MAIN_SCRIPT_CS_PATH
python genChromState.py $1
cd $MAIN_SCRIPT_CTCF_PATH
python genCTCFbinding.py $1

echo '>>>> The generated input files are located at `./runMolecularDynamics/inputFiles/epig_input/`'
echo
