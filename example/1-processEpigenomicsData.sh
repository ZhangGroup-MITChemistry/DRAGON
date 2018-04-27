#!/bin/bash

echo '''
****** Process Epigenomics Data ******'''

echo '
The chromatin is modeled as beads on a string. '
echo -n 'Each bead is assigned with a chromatin state, '
echo -n 'and will be labeled as a CTCF-binding site in a given orientation '
echo 'depending on the underlying Chip-Seq signal.
'
echo
echo '>>>> Chromatin states
'
echo -n 'ChromHMM is used to process epigenomics data of genome-wide histone modifications '
echo 'from six cell types and define the chromatin states. '
echo 'ChromHMM must be installed to run the script located at `./processEpigenomicsData/chromatinStates/run_chromHMM_5kb_6celltypes_15states.sh`. '
echo 'Detailed installation instructions for ChromHMM can be found on `http://compbio.mit.edu/ChromHMM/`. '
echo 'See `./processEpigenomicsData/README.md` for detailed information.'
echo
echo '>>>> NarrowPeak and motif files for CTCF-binding
'
echo 'The CTCF-binding sites are defined using the ChIP-Seq NarrowPeak binding profiles of CTCF together with Rad21 (Cohesin subunit). '
echo 'A near binding peak of Cohesin to the binding peak of CTCF is required to determine a CTCF-binding site. '
echo 'And the orientation of the CTCF-binding site is further determined using CTCF-binding motifs. '
echo 'See `./processEpigenomicsData/README.md` for more details on downloading and processing these files.'
echo
echo '>>>> The processed files are located in the individual folders in `./processEpigenomicsData/`.'
echo