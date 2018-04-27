#!/usr/bin/env python
import sys
sys.path.append('../src/md/')
from Ipt_module import *
from Params import *
Params()

from getSettings import getSettings
from CreateLAMMPSFile import *
from ProcessEpigData import *
from submitJobs import *

if __name__ == '__main__':

	# ---- input settings ---- #
	celltype,runnum,\
	nNode,ncpu,ptn,simtime,lmpsdir,\
	nearCtcfThreshold,runStep,\
	bind_flxb,cap,\
		chrom_lst=getSettings(sys.argv[1:])

	# ---- epigenomics data ready or not ---- #
	Ped=ProcessEpigData(celltype,chrom_lst,bind_flxb,cap)
	Ped.procChromState()
	Ped.procCTCFSites()

	# ---- cluster computing or locally ---- #
	clus_opt=raw_input('''  >> Computing clusters available?[y/n] ''')

	for chrId in chrom_lst:
		for runId in xrange(runnum):
			# ---- prepare simulation files ---- #
			Clf=CreateLAMMPSFile(celltype,chrId,runId,\
								nNode,ncpu,ptn,simtime,lmpsdir,\
									nearCtcfThreshold,runStep)
			Clf.createLAMMPSDataFile()
			Clf.createLAMMPSInputFile()

			if clus_opt == 'y':
				# ---- create and submit the job to the cluster ---- #
				Clf.createJobScript()
				submitJobs(celltype,chrId,runId); time.sleep(5);
			else:
				# ---- create the bash script to run locally ---- #
				Clf.createLocalBash()
	print
