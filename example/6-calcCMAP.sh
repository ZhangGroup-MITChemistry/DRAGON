#!/bin/bash

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
MAIN_SCRIPT_PATH=$SHELL_FOLDER/../analyzeChromatinConformation/contactMap/
MAIN_SCRIPT=$MAIN_SCRIPT_PATH/calContactMap.py
MAIN_SCRIPT_PATH_VIS='./analyzeChromatinConformation/contactMap/'

echo '''
****** Calculate the contact map ******'''
echo '
>>>> This script is to calculate the simulated contact map for chromosome 1, GM12878.'

PY_SCRIPT=""
LINE_ARR1=({1..21})
LINE_ARR2=({35..36})
LINE_ARR3=({38..41})
LINE_ARR=("${LINE_ARR1[@]}" "${LINE_ARR2[@]}" "${LINE_ARR3[@]}")
for NUM in ${LINE_ARR[@]};do
	LINE_RAW="$(sed -n "$NUM"p "$MAIN_SCRIPT")"
	if [ $NUM -eq 21 ];then
		LINE=$'	runnum=1\n'
	elif [ $NUM -eq 35 ];then
		LINE=${LINE_RAW:1}
	else
		LINE="$LINE_RAW"
	fi
	PY_SCRIPT="$PY_SCRIPT
$LINE"
done

cd $MAIN_SCRIPT_PATH
python -c "$PY_SCRIPT"

echo -n '>>>> The generated contact map is located at `./analyzeChromatinConformation/contactMap/cmap/`, '
echo "and please refer to MATLAB code: $MAIN_SCRIPT_PATH_VIS/visContactMap.m to visualize the generated contact map."
echo
