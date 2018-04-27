## Build orientation specific CTCF-binding sites

We use ChIP-Seq narrow-peak binding profiles for both CTCF and Rad21 (Cohesin subunit) to determine the list of loop forming CTCF-binding sites. For a CTCF site to be included in this list, there must be a Rad21 peak nearby as well. We encourage the interested reader to look into the supporting information of our manuscript and the python scripts provided below. The orientation of each CTCF-binding site is defined using CTCF-binding motifs as detailed below.


### NarrowPeak files

CTCF and Rad21 narrow-peak files can be downloaded by running the following command:
```
cd ./raw.narrowPeak/narrowPeak/; ./download.sh
```
These narrow-peak files are further processed by the python program [`./ctcfBindingSites/processingNarrowPeak.py`](./ctcfBindingSites/processingNarrowPeak.py) into txt files that contain only the genomic positions of individual binding sites for a given chromosome.

Usage:
```
python processingNarrowPeak.py [-C Celltype] [-t transcriptional_factor] [-c chromosome_id_array]
```
by default is calculating: Gm12878, chromosome 1.

**[Celltype]** can be selected from the following list:  
>Gm12878  
>H1hesc  
>Hela  
>Hepg2  
>Huvec  
>K562  

**[transcriptional_factor]** can be selected from the following list: 

> ctcf
> rad21  

**[chromosome_id_array]** can be any non-repeated subset selected from:  

>1 ~ 22

This script produces txt files that list the binding peaks for the sprecific cell type and chromosomes in the following format:

```
start_position(bp)     	end_position(bp)
114889141               114889474
150951988               150952271
225662661               225662969
145059036               145059036 
```
By default, the output of the binding narrow peak located in the folder `ctcfBindingSites/raw.narrowPeak/[Celltype]/ctcf(rad21)/`.

### Motif files

The motif files that are used to determine the orientation of CTCF-binding sites are from the following references:

>Rao, Suhas S.P. et al. *Cell* **159**, 1665-1680 (2014).  
>Kheradpour, P. & Kellis, M. *Nucleic Acids Res.* **42**, 2976-2987 (2014).

raw data for motifs can be downloaded by running the following command:

```
cd ./motif_file/motifs/; ./download.sh
```

The motif file for the first reference will not be downloaded by executing this command, but the link of the file is provided in [`./motif_file/motifs/download.sh`](./motif_file/motifs/download.sh).

The downloaded motif files is further processed by the python program [`./motif_file/prepareMotif.py`](./motif_file/prepareMotif.py) into txt files that contain only the genomic position and orientation of individual motifs for a given chromosome.

Usage: 
```
python prepareMotif.py [-m motif_file] [-p motif_folder_name_option] [-c chromosome_id_array]
```
**[motif_file]** is the name of the downloaded motif file.

**[motif_folder_name_option]** is the name of the processed folder.

**[chromosome_id_array]** can be any non-repeated subset selected from:  

> 1 ~ 22

This script produces txt files that list the positions and orientations of motifs for the specific chromosomes in the following format:

```
motif_position(bp)     	motif_orientation
12301               	-
13063               	-
15633               	-
17094               	+ 
```

By default is calculating: motif file from the first reference abovementioned and chromosome 1 (with the folder name as 'lbm'). The other two motif folders are calculated based on the motif files from the second reference abovementioned.
