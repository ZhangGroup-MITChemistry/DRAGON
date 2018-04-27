from Ipt_module import *
from Params import *
Params()

import lammps_tools as lmp


class CreateLAMMPSFile():

	_same_mol_flag_ideal = 0	
	_paramsFolder = '%s/../src/md/lmps_input/'%glb_path
	_lmpsTemplate = _paramsFolder+'lammps_template.in'

	def __init__(self,celltype,chrId,runId,nNode,ncpu,ptn,simtime,lmpsdir,nearCtcfThreshold,runStep):

		self.celltype           = celltype
		self.chrId              = chrId
		self.runId              = runId
		self.nNode              = nNode
		self.ncpu               = ncpu
		self.ptn                = ptn
		self.simtime            = simtime
		self.lmpsdir            = lmpsdir
		self.nearCtcfThreshold  = nearCtcfThreshold
		self.runStep            = runStep


	def assignType(self,seqfile,hp):
	# ----  assign the states based on seqfile 
		hp.atom_types = []
		seq = []
		for ics in range(1,ncs+1):
			hp.add_atom_type(1.0, None, None)
		for line in fi.input(seqfile):
			seq.append(int(line.split()[1]))
		for ia, atom in enumerate(hp.atoms):
			atom['atom_type_i'] = seq[ia]


	def createLAMMPSDataFile(self):
		global _gSta
		global _gEnd
		_gSta = chr_region[self.chrId-1,1]
		_gEnd = _gSta+25
		_ctcfFile 	= '%s/inputFiles/epig_input/ctcfSites/' \
			'%s/%s_chr%d_ctcf_position_From%dMbTo%dMb.txt'\
			%(glb_path,self.celltype,self.celltype,self.chrId,_gSta,_gEnd)
		_lmpsFolder	='%s/inputFiles/lmps_input/%s/'%(glb_path,self.celltype)

		# ----  specify bond and angle
		hp = lmp.Data()
		hp.add_bond_type(bond_coeffs,comment)
		hp.add_angle_type(angle_coeffs,comment)
		hp.read_from_file('%s/data.chromosome.init'%self._paramsFolder)
	 
		# ----  assign atom types
		self.assignType(_ctcfFile,hp)
		
		if not os.path.exists(_lmpsFolder):
			os.makedirs(_lmpsFolder)
		hp.write_to_file('%s/data.chromosome.chr%d'%(_lmpsFolder,self.chrId),ellipsoidFlag=0)		


	def createLAMMPSInputFile(self):
		global _rundir
		_rundir = "%s/run_folder/%s/chr%d/run%02d/"%(glb_path,self.celltype,self.chrId,self.runId)
		_csFile='%s/inputFiles/epig_input/chromStates/' \
			'%s/%s_chr%d_chromatin_states_From%dMbTo%dMb.txt'\
				%(glb_path,self.celltype,self.celltype,self.chrId,_gSta,_gEnd)
		_ctcfIndFile = '%s/inputFiles/epig_input/ctcfSites/' \
			'%s/%s_chr%d_ctcf_index_From%dMbTo%dMb.txt'\
			%(glb_path,self.celltype,self.celltype,self.chrId,_gSta,_gEnd)

		if not os.path.exists(_rundir):
			os.makedirs(_rundir)

		# ----  lammps input file with custom random seed
		_inFile = _rundir + "in.chromosome"
		in_tmp = fi.input(self._lmpsTemplate)
		pf = open(_inFile, 'w')
		pf.write('variable        rseed equal   %d\n'%(4928459+self.runId) )
		for line in in_tmp:
			if line[0:9] == 'read_data':
				pf.write(
'''read_data       ../../../../inputFiles/lmps_input/%s/data.chromosome.chr%d\n'''%(self.celltype,self.chrId))
			elif line[0:10] == 'pair_style':
				pf.write(
'''pair_style      hybrid/overlay table linear 10000 tanhlr/cut/ideala 6.0 %d 15 %s/ucs_chrom.txt %s tanhlr/cut/ideal 6.0 %s/uctcf_chrom.txt %s %d\n'''\
%(self._same_mol_flag_ideal,\
self._paramsFolder,_csFile,\
self._paramsFolder,_ctcfIndFile,self.nearCtcfThreshold))
			elif line[0:16] == 'pair_coeff_softc':
				pf.write(
'''pair_coeff       * * table %s/soft_core_lj_4kT.table soft_core_lj 1.12\n'''%(self._paramsFolder))
			elif line[0:3] == 'run':
				pf.write('run             %d'%(self.runStep))
			else:
				pf.write(line)

		pf.close()


	def createJobScript(self):
		# ----  create pbs job script (parallel run)
		pbsFile = _rundir + "job.pbs"
		pf = open(pbsFile, 'w')
		pbs = '''#!/bin/bash

#SBATCH --job-name=%s_c%d
#SBATCH -N %d
#SBATCH -n %d
#SBATCH --partition=%s
#SBATCH --no-requeue
#SBATCH --time=%d:00:00
#SBATCH --export=ALL

module load gcc
module add mvapich2/gcc

lammpsdir="%s"
mpirun -np %d $lammpsdir/lmp_openmpi -in in.chromosome
'''%(self.celltype,self.chrId,self.nNode,self.ncpu,self.ptn,self.simtime,self.lmpsdir,self.ncpu)
		pf.write(pbs)
		pf.close()
		print('''   > Job for %s, chromosome %d, parrallel running %02d is processed for submission.'''\
									%(self.celltype,self.chrId,self.runId))

	def createLocalBash(self):
		# ----  create local bash script (serial run)
		localtime = time.asctime( time.localtime(time.time()) )
		runFile = _rundir+"run.sh"
		pf = open(runFile, 'w')
		pbs = '''# ---- Chromosome simulation locally in serial ---- #
# ---- Created Time: %s ---- #

#!/bin/bash
%s/lmp_openmpi -in in.chromosome'''%(localtime,self.lmpsdir)
		pf.write(pbs)
		pf.close()

		_cmd = 'chmod 744 %s'%(runFile)
		q = Popen(_cmd, shell=True, stdout=PIPE)
		q.communicate()
		print('''   > Local simulation bash script is generated located at: %s/run_folder/%s/chr%d/run%02d/run.sh'''\
									%(glb_path,self.celltype,self.chrId,self.runId))
