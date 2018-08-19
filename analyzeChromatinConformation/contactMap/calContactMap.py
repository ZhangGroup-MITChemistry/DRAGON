# build up parallel simulations for calculation of cmaps
import sys
sys.path.append('../../src/cmap/PYTHON/')
from Ipt_module import *
from Params import *
Params()

from getSettings import getSettings
from processingJobScript import *
from calMapLocal import calMapLocal
from checkStatus import checkStatus
from combineMaps import combineMaps


if __name__ == '__main__':

	# ---- input settings ---- #
	celltype,runnum,\
	jobname,usrname,ptn,\
		chrom_lst=getSettings(sys.argv[1:])

	# ---- cluster computing or locally ---- #
	clus_opt=raw_input("   > Computing clusters available?[y/n] ")
	
	# ---- cluster available ---- #
	if clus_opt == 'y':
		# ---- prepare the job script for calculating cmaps ---- #
		processingJobScript(celltype,runnum,jobname,ptn,chrom_lst)

		# ---- check the job status ---- #
		checkStatus(usrname,jobname)

	# ---- compute locally ---- #
	else:
		calMapLocal(celltype,runnum,chrom_lst)

	# ---- combine the parallel cmaps to ensemble average ---- #
	for chrId in chrom_lst:
		combineMaps(celltype,runnum,chrId)
	print
