## LAMMPS Files for DRAGON 

This folder contains our custom modifications to the LAMMPS software package to enable structure prediction for chromatin. These files are compiled together with the main package using the script [`LAMMPS.sh`](../LAMMPS.sh). See the Installation section of the main [README](../../README.md) page for details.

### Description

[`pair_tanhlr_cut_ideala.cpp`](./pair_tanhlr_cut_ideal.cpp) is a LAMMPS source code for the sequence dependent chromatin state specific interaction between pairs of loci. This function takes the [sequence of chromatin states](../../runMolecularDynamics/inputFiles/epig_input/chromStates/) and the [optimized contact energies for chromatin states interactions](../md/lmps_input/ucs_chrom.txt) as input. 

[`pair_tanhlr_cut_ideala.h`](./pair_tanhlr_cut_ideal.h) is the header file for [`pair_tanhlr_cut_ideala.cpp`](./pair_tanhlr_cut_ideal.cpp).

[`pair_tanhlr_cut_ideal.cpp`](./pair_tanhlr_cut_ideal.cpp) is a LAMMPS source code for interactions between pairs of genomic loci enclosed by convergent CTCF pairs. This function takes the [sequence of CTCF indices](../../runMolecularDynamics/inputFiles/epig_input/ctcfSites/), the [optimized contact energies for CTCF interactions](../md/lmps_input/uctcf_chrom.txt) and the threshold for nearest CTCF pair that mimics the finite processivity in the extrusion model as input. 

[`pair_tanhlr_cut_ideal.h`](./pair_tanhlr_cut_ideal.h) is the header file for [`pair_tanhlr_cut_ideal.cpp`](./pair_tanhlr_cut_ideal.cpp).

For detailed descriptions of the energy functions associated with above source code, please refer to the [Supplementary Material](https://www.biorxiv.org/highwire/filestream/86852/field_highwire_adjunct_files/0/282095-1.pdf) of the [manuscript](https://www.biorxiv.org/content/early/2018/03/15/282095).

[`force.cpp`](./force.cpp) is an ordinary LAMMPS source code. It is put here simply to avoid the possible incompatibility between different versions of LAMMPS.

[`Makefile.openmpi`](./Makefile.openmpi) is a make file that is used to compile the parallel version of LAMMPS.
