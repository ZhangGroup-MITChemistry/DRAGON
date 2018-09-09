## Chromatin Structure Prediction using Molecular Dynamics Simulation 

The python script [`main.py`](./main.py) provides a streamlined interface for predicting chromosome structures. It combines the three steps listed in the Usage section of the main [README](../README.md) page, to extract epigenomics input, to build LAMMPS input, and to run simulations. As detailed below, using this script, one can take advantage of parallel computing to simulate chromosome structures from a variety of cell types.

Usage:
```
python main.py --Cell <Celltype> --chrom <chromosome id> --region <chromosome region> [--lmpsdir <Lammps dir>] [--step <simulation steps>] [--njob number_of_jobs] [--nNode <number of Nodes per job>] [--ncpu <number of cpu per job>] [--ptn <partition>] [--time <simulation time>]  
```

The items enclosed in the brackets [] are optional parameters that are used to produce scripts required by the [Slurm job scheduling system](https://slurm.schedmd.com/) to perform simulations on a computer cluster. 

**[--Cell]** specifies the cell type of interest. The list of supported cell types are (case sensitive):

>Gm12878  
>H1hesc  
>Hela  
>Hepg2  
>Huvec  
>K562

**[--chrom]** specifies the list of chromosomes to be simulated. The accepted format is:
> 1,10,19,21

Any number from 1 to 22 can be included in this list.

**[--region]** specifies the genomic regions of interest for each chromosome specified in [--chrom]. The accepted format is:

> [20,45],[88,113],[34,59],[20,45]

The number of entries included in the list [] must be same as those in [--chrom]. If only one region is specified, it will be applied to all the chromosomes specified in [--chrom]. The default is the genomic region from 20Mb to 45Mb. The unit of these numbers is million base pairs. We recommend simulating chromosome regions of 25Mb or less for computational efficiency and accuracy, but simulations for longer regions are also feasible. 

**[--lmpsdir]** specifies the path to the compiled LAMMPS binary file, and the accepted format is:

> `/path-to-LAMMPS-folder/src/`. 

See [README](../README.md) for more details on installing LAMMPS. If the LAMMPS is successfully installed, the default path to be directed to the compiled binary file where the LAMMPS is installed.

**[--step]** specifies the number of MD steps that will be performed. The default value is set to be 40 million steps.  

**[--njob]** specifies the number of independent MD simulations that will be performed for each chromosome.  By default, it is set to be 1.

**[--nNode]** **[--ncpu]** **[--ptn]** **[--time]** should be chosen based on the configuration of the computer cluster. By default, they are set to be 1 Node, 14 cpus, 'mit', and 48hrs.  

The complete list of options is also available via the command:
```
python main.py -h
```
