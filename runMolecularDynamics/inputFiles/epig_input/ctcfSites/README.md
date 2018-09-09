## Extracting CTCF-binding sites for the chromatin region of interest
This script produces orientation specific CTCF-binding sites for the 25Mb long chromatin regions defined in the [txt file](../../../../src/chr_region.txt).  

### Output
Two txt files will be produced from this script.

#### 1) list of orientation specific CTCF binding sites
This file contains two columns, with the first one correponding to the polymer bead ID and the second for type of CTCF-binding sites. Four CTCF types are defined.

>1 for 5-3 oriented binding sites  
>2 for 3-5 oriented binding sites  
>3 for non-binding sites  
>4 for multiple binding sites with both orientations  

#### 2) index file

This file contains two columns that provide the location of nearest 5-3 and 3-5 CTCF binding sites for a given polymer bead respectively.

### Usage:

```
python genCTCFbinding.py --Cell <Celltype> --bindflex <motif_match_flexibility> --cap <CTCF-cohesin_nearest_dist> --chrom <chromosome_id_array>
```
or default:
```
python genCTCFbinding.py
```
Note items in [] are optional. By default is calculating: Gm12878, chromosome 1.  

**[--Cell]** can be selected any one from the following list (case sensitive):

>Gm12878  
>H1hesc  
>Hela  
>Hepg2  
>Huvec  
>K562

**[--bindflex]** is the parameter that indicates the region for the motif matching. The default value is set to be 100bp, which means that if a motif is found within 100bp of the binding site, the orientation of the binding site was then assigned based on the DNA sequence of that motif.  

**[--cap]** is the parameter that indicates the requirement for the nearest genomic distance of CTCF with Rad21. The default value is set to be 50bp.  

**[--chrom]** can be any non-repeated subset selected from 1 to 22. The accepted format is:

> 1 10 19 21

The CTCF-binding sites input to the model is processed based on the NarrowPeak binding profiles and motifs located in this [folder](../../../../processEpigenomicsData/ctcfBindingSites/). The generated files as the input of modeling are located in the folder: `./[Celltype]/`.

The manual would be available by executing:

```
python genCTCFbinding.py -h
```
