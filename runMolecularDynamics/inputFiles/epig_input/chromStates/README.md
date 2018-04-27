## Extract chromatin states for the chromatin region of interest

This script processes the [output](../../../../processEpigenomicsData/chromatinStates/OUTPUTSAMPLE_5kb_6celltype_15states/) from ChromHMM to extract chromatin states for the 25Mb long chromatin regions defined in the [txt file](../../../../src/chr_region.txt).  

### Output

A two column txt file will be produced, with the first column indicates the polymer bead ID and the second for chromatin state ID. 

### Usage:

```
python genChromState.py [-C Celltype] [-c chromosome_id_array]
```
or default:
```
python genChromState.py
```
Note items in [] are optional. By default is calculating: Gm12878, chromosome 1.

**[Celltype]** can be selected from the following list (case sensitive):
>Gm12878  
>H1hesc  
>Hela  
>Hepg2  
>Huvec  
>K562

**[chromosome_id_array]** can be any non-repeated subset selected from:
>1 ~ 22

The manual would be available by executing:
```
python genChromState.py -h
```
