## Chromatin Structure Prediction using Molecular Dynamics Simulation 

The python script [`main.py`](./main.py) provides a simple interface to perform molecular dynamics simulations. It combines together the three steps listed in the Usage section of the main [README](../README.md) page, to extract epigenomics input, to build LAMMPS input, and to run simulations. 

Furthermore, this script can set up multiple parallel simulations for improved sampling and computational efficiency. 

Usage:
```
python main.py [-C Celltype] [-n runnumber] [-N number_of_Node] [-p number_of_cpu] [-i partition] [-t simulation_time] [-l Lammps_dir] [-d near_CTCF_threshold] [-r simulation_steps] [-b motif_match_flexibility] [-a CTCF-cohesin_nearest_dist] [-c chromosome_id_array]
```
or default:
```
python main.py
```
Note items in [] are optional. By default is calculating: Gm12878, chromosome 1, 8 parallel running.  

**[Celltype]** can be selected from the following list (case sensitive):
>Gm12878  
>H1hesc  
>Hela  
>Hepg2  
>Huvec  
>K562

**[runnumber]** specifies the number of parallel running. By default, the value is 8.  

**[number_of_Node]** **[number_of_cpu]** **[partition]** **[simulation_time]** should be specified based on the computational capability. By defaul, the values are set to be 1 Node, 14 cpus, 48hrs at maximum.  

**[Lammps_dir]** should be specified as `/path-to-LAMMPS-folder/src/`. The default path has already been directed to the installed LAMMPS. If the LAMMPS package has not been installed yet, see [README](../README.md) for more details.  

**[near_CTCF_threshold]** is the parameter that indicates the number of convergent CTCF pairs that are separated by no more than a finite number of CTCF-binding sites in both orientations to mimic the finite processivity of cohesin molecules. The default value is set to be 4.  

**[simulation_steps]** indicates the number of MD steps. The default value is set to be 40 million steps.  

**[motif_match_flexibility]** is the parameter that will be used during the processing of CTCF-binding sites. This parameter indicates the region for the motif matching. See [README](./inputFiles/epig_input/ctcfSites/README.md) for detail.  

**[CTCF-cohesin nearest dist]** is the parameter that will be used during the processing of CTCF-binding sites. This parameter indicates the requirement for the nearest genomic distance of CTCF with Rad21. See [README](./inputFiles/epig_input/ctcfSites/README.md) for detail.  

**[chromosome_id_array]** can be any non-repeated subset selected from:
>1 ~ 22

The manual would be available by executing:
```
python main.py -h
```
