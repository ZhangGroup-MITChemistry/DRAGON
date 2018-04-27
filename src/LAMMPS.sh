#!/bin/bash

git clone -b stable https://github.com/lammps/lammps.git lammps

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
LMPS_FOLDER_SRC=$SHELL_FOLDER/../lammps/src/
LMPS_SRC=$SHELL_FOLDER/lammps/

cp $LMPS_SRC/Makefile.openmpi $LMPS_FOLDER_SRC/MAKE/
cp $LMPS_SRC/force.cpp $LMPS_FOLDER_SRC/
ln -s $LMPS_SRC/pair_tanhlr_cut_ideal.cpp $LMPS_FOLDER_SRC/pair_tanhlr_cut_ideal.cpp
ln -s $LMPS_SRC/pair_tanhlr_cut_ideal.h $LMPS_FOLDER_SRC/pair_tanhlr_cut_ideal.h
ln -s $LMPS_SRC/pair_tanhlr_cut_ideala.cpp $LMPS_FOLDER_SRC/pair_tanhlr_cut_ideala.cpp
ln -s $LMPS_SRC/pair_tanhlr_cut_ideala.h $LMPS_FOLDER_SRC/pair_tanhlr_cut_ideala.h

echo '''
>>>> Connected with the LAMMPS src files.'''
echo '''
>>>> Starting to compile LAMMPS ......
'''
cd $LMPS_FOLDER_SRC
make clean-all
make yes-molecule
make yes-class2
make -j 4 openmpi
