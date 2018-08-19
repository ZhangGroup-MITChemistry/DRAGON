#!/bin/bash

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
MAIN_SCRIPT_PATH=$SHELL_FOLDER/../runMolecularDynamics/
MAIN_SCRIPT=$MAIN_SCRIPT_PATH/main.py

echo '''
****** Build LAMMPS Input ******'''
echo '
>>>> This script is to build the LAMMPS input scripts for chromosome 1, GM12878.'

PY_SCRIPT=""
LINE_ARR1=({1..21})
LINE_ARR2=({30..37})
LINE_ARR3=({45..47})
LINE_ARR=("${LINE_ARR1[@]}" "${LINE_ARR2[@]}" "${LINE_ARR3[@]}")
for NUM in ${LINE_ARR[@]};do
	LINE_RAW="$(sed -n "$NUM"p "$MAIN_SCRIPT")"
	if [ $NUM -eq 21 ];then
		LINE=$'	runnum=1\n'
	elif [ $NUM -eq 45 ];then
		LINE=${LINE_RAW:1}
	else
		LINE="$LINE_RAW"
	fi
	PY_SCRIPT="$PY_SCRIPT
$LINE"
done

cd $MAIN_SCRIPT_PATH
python -c "$PY_SCRIPT"

echo '>>>> The generated input files are located at `./runMolecularDynamics/run_folder/Gm12878/chr1/run00/`'
echo