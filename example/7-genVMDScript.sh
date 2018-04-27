#!/bin/bash

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
MAIN_SCRIPT_PATH=$SHELL_FOLDER/../analyzeChromatinConformation/visStructure/

echo '''
****** Generate VMD Script ******'''
echo '''
>>>> This script is to generate the VMD script for the visualization of 3D structure for chromosome 1, GM12878.'''

cd $MAIN_SCRIPT_PATH
python VMDmain.py

echo -n '>>>> The generated script is located at `./analyzeChromatinConformation/visStructure/`, '
echo 'and please refer to the README file located in that folder for detailed instruction of visualization steps with VMD.'
echo
