# ![DRAGON logo](https://github.com/qiyf/images/blob/master/logo2.png)

[![GitHub shield](https://img.shields.io/badge/GitHub-DRAGON-orange.svg?style=flat)](https://github.com/ZhangGroup-MITChemistry/DRAGON) [![](https://img.shields.io/badge/version-latest-yellow.svg)](https://github.com/ZhangGroup-MITChemistry/DRAGON) [![bioRxiv shield](https://img.shields.io/badge/bioRxiv-1709.01233-green.svg?style=flat)](https://www.biorxiv.org/content/early/2018/03/15/282095) [![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Contents

- [Overview](#overview)

- [System Requirements](#system-requirements)

- [Installation](#installation)

- [Usage](#usage)

  + [Process Epigenomics Data](#i-process-epigenomics-data)

  + [Run Molecular Dynamics Simulation](#ii-run-molecular-dynamics-simulation)

  + [Analyze Chromatin Conformation](#iii-analyze-chromatin-conformation)

- [License](./LICENSE)

## Overview

DRAGON is a software package to enable De novo, and RAtional prediction of Genome organizatiON. It provides an implementation of the model proposed in the [manuscript](https://www.biorxiv.org/content/early/2018/03/15/282095) to simulate chromatin structure and dynamics. With DRAGON, one can predict the structure of a 25Mb long chromatin region from a variety of cell types using genome-wide profiles of histone modifications and CTCF molecules. 

The package is mainly written in [Python](https://www.python.org/), and it streamlines all the necessary steps to process [epigenomics data](./processEpigenomicsData/), to perform [molecular dynamics simulations](./runMolecularDynamics/) and to analyze [predicted conformational ensemble](./analyzeChromatinConformation/) for the chromatin. 

## System Requirements

### Hardware

DRAGON requires a single cpu (serial) or multiple cpus (parallel) with standard RAM to run molecular dynamics simulations. For a decent level of preformance, we recommend computing node with 14+ cores, 2.6+ GHz/core. All simulations in our [manuscript](https://www.biorxiv.org/content/early/2018/03/15/282095) are performed on an Intel Xeon E5-2690 v4 2.6Ghz node with 14 cores. 

### Operating System

The package is supported for Linux operating systems. It had been tested on the system Linux CentOS 6.8.

## Installation

DRAGON can be installed by running the following command:
```
git clone https://github.com/ZhangGroup-MITChemistry/DRAGON.git
```
or by downloading the zip file with the link:

[https://github.com/ZhangGroup-MITChemistry/DRAGON/archive/master.zip](https://github.com/ZhangGroup-MITChemistry/DRAGON/archive/master.zip)  

DRAGON uses [LAMMPS](http://lammps.sandia.gov/), a molecular dynamics software package, for simulating chromatin structures. LAMMPS, together with our custom modifications, can be compiled and installed with the following command:

```
./src/LAMMPS.sh
```

Note that the [GCC](https://gcc.gnu.org/) compiler needs to be installed beforehand and an environment of [OpenMPI](https://www.open-mpi.org/) is needed to compile the parallel version of LAMMPS. 

On a "normal" desktop computer, the typical install time including download and compile steps is ~ 5 minutes. 

## Usage

DRAGON models the chromatin as a coarse-grained bead-spring polymer, with each bead corresponding to a five kb genomic segment.  This coarse-grained polymer is made cell type and chromosome specific by assigning each bead with a chromatin state. The polymer bead will also be labeled as an orientation dependent CTCF binding site if there is a strong binding signal in the corresponding region. With its underlying biochemistry specified, the structure of the chromatin can be predicted by simulating the sequence-specific potential energy function parameterized in our [manuscript](https://www.biorxiv.org/content/early/2018/03/15/282095) using LAMMPS. See the flowchart below for an illustration of the different steps for chromatin structure prediction.

![Flow chart](https://github.com/qiyf/images/blob/master/flow_chart.png)

We further provide step-by-step instructions below to simulate the structure of chromosome 1 from GM12878 cells. All the executable scripts are provided in the [`./example/`](./example/) folder. 

A [python script](./runMolecularDynamics/main.py) is also provided to set up simulations for other chromosomes from different cell types. It streamlines all the steps below together and can launch multiple parallel simulations. See [README](./runMolecularDynamics/README.md) for its detailed usage.


### I) Process Epigenomics Data

Before starting any structure predictions, we need to learn the chromatin states from genome-wide histone modification profiles and identify the genomic location and orientation of CTCF binding sites. 

```
./example/1-processEpigenomicsData.sh
```

This script provides detailed instructions on how to process epigenomics data using [ChromHMM](http://compbio.mit.edu/ChromHMM/) and custom python scripts. 

### II) Run Molecular Dynamics Simulation

To start simulating chromatin structures, the following steps are necessary in order to incorporate the processed epigenomics input into data formats recognized by LAMMPS.


#### Select a 25Mb chromatin region

First, one needs to select a 25Mb long chromatin region of interest (the default setting for this example running is chr1:20-45Mb from GM12878 cells) by running the following script 

```
./example/2-selectChromatinRegion.sh
```

DRAGON currently can only simulate chromatin regions with a fixed length of 25Mb, but generalization to whole chromosomes is straightforward. 

This script produces a [txt file](./src/chr_region.txt) that lists the region of interested for each chromosome in the following format:
```
chromosome_id     start_position(Mb)      end_position(Mb)  
1                 20                      45  
2                 20                      45  
3                 20                      45  
4                 20                      45   
```

If a different chromatin region is desired, one can either directly modify the generated chromatin region [txt file](./src/chr_region.txt) or modified the parameters in the original script [`./example/2-selectChromatinRegion.sh`](./example/2-selectChromatinRegion.sh) and then regenerate the file.

#### Extract epigenomics input

Second, one needs to extract chromatin states and CTCF binding sites for the selected chromatin region from results produced in **section I)**.

```
./example/3-extractEpigenomicsInput.sh
```

Three txt files are generated to provide the [chromatin state of each polymer bead](./runMolecularDynamics/inputFiles/epig_input/chromStates/Gm12878/Gm12878_chr1_chromatin_states_From20MbTo45Mb.txt), the [CTCF binding potency of each polymer bead](./runMolecularDynamics/inputFiles/epig_input/ctcfSites/Gm12878/Gm12878_chr1_ctcf_position_From20MbTo45Mb.txt), and the [location of nearest CTCF binding sites for each polymer bead](./runMolecularDynamics/inputFiles/epig_input/ctcfSites/Gm12878/Gm12878_chr1_ctcf_index_From20MbTo45Mb.txt). 
See [Chromatin States README](./runMolecularDynamics/inputFiles/epig_input/chromStates/README.md) and [CTCF-binding Sites README](./runMolecularDynamics/inputFiles/epig_input/ctcfSites/README.md) for details on file formats and data extraction.

#### Build LAMMPS input

Third, one needs to incorporate the epigenomic inputs produced from the step above into file formats recognized by LAMMPS for molecular dynamics simulation.

```
./example/4-buildLammpsInput.sh
```

This script produces a [topology file](./runMolecularDynamics/inputFiles/lmps_input/Gm12878/data.chromosome.chr1) that stores the Cartesian coordinates of each polymer beads and the connectivity among polymer beads, an [input file](./runMolecularDynamics/run_folder/Gm12878/chr1/run00/in.chromosome) that instructs the specifics of the molecular dynamics simulation, and a [bash script](./runMolecularDynamics/run_folder/Gm12878/chr1/run00/run.sh) to execute LAMMPS. 

#### Run simulation

With all the preparation steps above, we are finally ready to predict chromatin structures by running the following script:

```
./example/5-runMD.sh
```

Simulated chromatin structures will be stored in a binary file (DUMP_FILE.dcd) in the folder [`./runMolecularDynamics/run_folder/Gm12878/chr1/run00/`](./runMolecularDynamics/run_folder/Gm12878/chr1/run00/). See below on how to use VMD to read this binary file and visualize chromatin structures. Note that this script only runs one molecular dynamics simulation using one CPU. We typically perform multiple simulations to improve conformational sampling and conduct each simulation with multiple CPUs to reduce simulation time. See the [python program](./runMolecularDynamics/main.py) on how to setup multiple parallel simulations. 


### III) Analyze Chromatin Conformation

We use [MATLAB](https://www.mathworks.com/products/matlab.html) to analyze contact maps and [VMD](http://www.ks.uiuc.edu/Research/vmd/) to visualize chromatin structures. Installation of these two software packages is highly recommended. See [contact map README](./analyzeChromatinConformation/contactMap/README.md) and [structure visualization README](./analyzeChromatinConformation/visStructure/README.md) for detailed instructions on usage. 

#### Calculate and visualize the contact map

To quantitatively compare predicted chromatin structures with genome-wide chromosome conformation capture experiments (Hi-C), one needs to first calculate the contact probability map by running the script:

```
./example/6-calcCMAP.sh
```

The [core program](./src/cmap/FORTRAN/cmap.f90) to calculate the contact map from trajectory files is written in FORTRAN and has been precompiled with ifort compiler. One can also modify the [script](./src/cmap/FORTRAN/compile.sh) to compile it with gfortan if ifort is not available. 

The calculated contact map is located at [`./analyzeChromatinConformation/contactMap/cmap/`](./analyzeChromatinConformation/contactMap/cmap/). To visualize the contact map, use the provided [MATLAB script](./analyzeChromatinConformation/contactMap/visContactMap.m). See [README](./analyzeChromatinConformation/contactMap/README.md) for more detailed instructions. 


#### Visualize the 3D structure

To render the predicted chromatin structures in 3D, run the following script:

```
./example/7-genVMDScript.sh
```

The script produces a [VMD file](./analyzeChromatinConformation/visStructure/vmdScript/VMDColor_Gm12878_chr1.vmd) written in TCL language that can be loaded into the software VMD to visualize chromatin structures. See [README](./analyzeChromatinConformation/visStructure/README.md) for detailed instruction of visualization steps with VMD.
