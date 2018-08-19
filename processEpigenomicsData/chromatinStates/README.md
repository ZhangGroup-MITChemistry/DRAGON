## Build Chromatin States with ChromHMM

We learned 15 chromatin states using 12 key histone modifications 
>H3K4me1  
>H3K4me2  
>H3K4me3  
>H3K9ac  
>H3K9me3  
>H3K27ac  
>H3K27me3  
>H3K36me3  
>Dnase  
>H2A.Z  
>H3K79me2  
>H4K20me1  

for the following six cell types :
>Gm12878  
>H1hesc  
>Hela  
>Hepg2  
>Huvec  
>K562 

Genome-wide histone modification profiles can be downloaded from [ENCODE](https://www.encodeproject.org) using the links provided in the [Extend Data Sheet](https://www.biorxiv.org/highwire/filestream/86852/field_highwire_adjunct_files/1/282095-2.xlsx) of the [manuscript](https://www.biorxiv.org/content/early/2018/03/15/282095).

The software [ChromHMM](http://compbio.mit.edu/ChromHMM/) is used to learn chromatin states from these genome-wide data and must be installed for the analysis below.

To learn chromatin states, run the following script:
```
./run_chromHMM_5kb_6celltypes_15states.sh
```
For the script to run successfully, one needs to modify 'path_to_ChromHMM' to the directory where the ChromHMM software is located, and 'path_to_bedfiles' to the directory where the downloaded histone modification files are located. 

Example outputs are provided in the folder [`./OUTPUTSAMPLE_5kb_6celltype_15states/`](./OUTPUTSAMPLE_5kb_6celltype_15states/) and [`./binarizedData_5kb_6celltype_15states/`](./binarizedData_5kb_6celltype_15states/). 

Note that if new cell types are added, be aware of corresponding the new chromatin states with the present-defined ones through the generated emission pattern. 
