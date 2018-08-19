# generate the job script and submit the job
from Ipt_module import *
from Params import *
Params()

def processingJobScript(celltype,runnum,jobname,ptn,chrom_lst):
# ---- generate the job script for calculating the contact map ---- #
	src_path = '%s/../../src/cmap/FORTRAN/'%glb_path
	ipt_path = '../../runMolecularDynamics/run_folder/'

	for chrId in chrom_lst:
		gSta = chr_region[chrId-1,1]
		gEnd = chr_region[chrId-1,2]
		sepDist = gEnd-gSta
		nbead = int(sepDist*Mb/resolution)
		# ---- default: calculate cmap of 80% region in the middle ---- #
		startb = nbead*0.08;
		endb = nbead*0.88;
		for runid in xrange(runnum):
			dcd_path  = '%s/%s/%s/chr%d/run%02d/'\
						%(glb_path,ipt_path,celltype,chrId,runid)
			cmap_path = './%s/chr%d/run%02d/'\
						%(celltype,chrId,runid)

			if not os.path.exists(cmap_path):
				os.makedirs(cmap_path)

			fo = open('%s/job_cg.pbs'%(cmap_path),'w')
			fo.writelines('''#!/bin/bash

#SBATCH --job-name=%s
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=%s
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --export=ALL

%s/cmap %s/DUMP_FILE.dcd %d %d %d %d'''\
%(jobname,ptn,src_path,dcd_path,cg_fac,startb,endb,startfr))
			fo.close()

			cmd = 'cd %s;sbatch job_cg.pbs;'%(cmap_path)
			q = Popen(cmd, shell=True, stdout=PIPE)
			q.communicate()
		print('''   > Job for calculating contact map for %s, chromosome %d is submitted.'''\
				%(celltype,chrId))
	print