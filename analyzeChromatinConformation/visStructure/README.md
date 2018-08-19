## Visualize the 3D structure

### Generate VMD script

To visualize the predicted chromatin structures, one first needs to run the python script [`VMDmain.py`](./VMDmain.py) to generate a PDB file, a PSF file, and a VMD script. 

The PDB file provides the Cartesian coordinates of the structure. The PSF file defines the topology of the chromatin and the connectivity between polymer beads. The VMD script instructs VMD to read these two files and render a 3D structure in which the chromatin is colored by chromatin states.

Usage:
```
python VMDmain.py [-C Celltype] [-c chromosome_id_array]
```
The manual would be available by executing:  
```
python VMDmain.py -h
```

### Visualize 3D structure with VMD

To load the files generated above, open the VMD software (available at [http://www.ks.uiuc.edu/Research/vmd/](http://www.ks.uiuc.edu/Research/vmd/)) and click under `Extension/TK Console` option. This will open up a TCL command window. 

First, move to the folder: [`./VMDScript/`](./VMDScript/) by typing the following in the TCL command window:
```
cd ./VMDScript/
```
Next, load the VMD script by typing:
```
play VMDColor_Gm12878_chr1.vmd
```
Now, a chromosome structure colored by chromatin states should be appearing in the main display window. 

To visualize the ensemble of predicted chromatin structures that is stored in the DUMP_FILE.dcd file, first type the following command in the TCL command window: 
```
cd ../../runMolecularDynamics/run_folder/Gm12878/chr1/run00/
```
This will move the directory to the folder where the simulation for chr1 from GM12878 was carried out. Change the path accordingly for other chromosomes. 

Then type the following command to load the dcd file:
```
animate read dcd DUMP_FILE.dcd 
```
where the trajectory would be loaded and appearing in the display window.
