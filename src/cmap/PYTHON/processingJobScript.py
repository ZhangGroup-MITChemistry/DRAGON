# generate the job script and submit the job
from Ipt_module import *
from Params import *
Params()

def processingJobScript(celltype,runnum,jobname,ptn,chrom_lst):
# ---- generate the job script for calculating the contact map ---- #
	for chrId in chrom_lst:
		for runid in xrange(runnum):
			dcd_path  = '%s/../../runMolecularDynamics/' \
				'/run_folder/%s/chr%d/run%02d/'\
				%(glb_path,celltype,chrId,runid)
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
#SBATCH --mem-per-cpu=2G
#SBATCH --time=4:00:00
#SBATCH --export=ALL

%s/../../src/cmap/FORTRAN/cmap %s/DUMP_FILE.dcd %d %d %d %d'''\
%(jobname,ptn,glb_path,dcd_path,cg_fac,startb,endb,startfr))
			fo.close()

			cmd = 'cd %s;sbatch job_cg.pbs;'%(cmap_path)
			q = Popen(cmd, shell=True, stdout=PIPE)
			q.communicate()
		print('''   > Job for calculating contact map for %s, chromosome %d is submitted.'''\
				%(celltype,chrId))
	print