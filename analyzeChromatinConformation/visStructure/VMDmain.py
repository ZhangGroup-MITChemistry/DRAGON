# generate VMD scripts
import sys
sys.path.append('../../src/vmd/')

from getSettings import getSettings
from VMD import *


if __name__ == '__main__':

	# ---- input settings ---- #
	celltype,chrom_lst = getSettings(sys.argv[1:])

	Vmd=VMD()
	for chrId in chrom_lst:
		
		# ---- generate individual pdb/psf files ---- #
		Vmd.genPdb(celltype,chrId)
		Vmd.genPsf(celltype,chrId)

		# ---- generate individual VMD scripts ---- #
		Vmd.genVMDScript(celltype,chrId)
	print
