from Ipt_module import *
from Params import *
Params()

def submitJobs(celltype,chrId,runId):

	# ----	iteratively submit the jobs ----
	rundir = "%s/run_folder/%s/chr%d/run%02d/"%(glb_path,celltype,chrId,runId)
	cmd = "cd %s; sbatch job.pbs;"%rundir
	q = Popen(cmd, shell=True, stdout=PIPE)
	q.communicate()