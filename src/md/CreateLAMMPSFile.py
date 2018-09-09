from Ipt_module import *
from Params import *

import lammps_tools as lmp
import polymer_tools as poly

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
		global _sepDist
		try:
			_gSta = chr_region[str(self.chrId)][0]
			_gEnd = chr_region[str(self.chrId)][1]
		except KeyError:
			with open('%s/../src/chr_region.txt'%glb_path,'r') as f:
				chr_reg = json.load(f)
				_gSta = chr_reg[str(self.chrId)][0]
				_gEnd = chr_reg[str(self.chrId)][1]
		_sepDist = _gEnd-_gSta

		global _ctcfFile
		_ctcfFile 	= 	'%s/inputFiles/epig_input/ctcfSites/' \
						'%s/%s_chr%d_ctcf_position_From%dMbTo%dMb.txt'\
						%(glb_path,self.celltype,self.celltype,self.chrId,_gSta,_gEnd)
		_lmpsFolder	=	'%s/inputFiles/lmps_input/%s/'%(glb_path,self.celltype)
		_dataFile	=	'%s/data.chromosome.init%d'%(self._paramsFolder,_sepDist)

		
		# ----	generate initial configuration through rondom walk
		if not os.path.exists(_dataFile):
			geninit = poly.GenInitConfig(_sepDist,self.lmpsdir,_ctcfFile)
			geninit.genInitConfig()

		# ----  specify bond and angle
		hp = lmp.Data()
		hp.add_bond_type(bond_coeffs,comment)
		hp.add_angle_type(angle_coeffs,comment)

		# ----  assign atom types
		hp.read_from_file(_dataFile)
		self.assignType(_ctcfFile,hp)

		lmp.make_folder(_lmpsFolder)
		hp.write_to_file('%s/data.chromosome.chr%d'%(_lmpsFolder,self.chrId),ellipsoidFlag=0)		


	def createLAMMPSInputFile(self):
		global _rundir
		_rundir = 	"%s/run_folder/%s/chr%d/run%02d/"%(glb_path,self.celltype,self.chrId,self.runId)
		lmp.make_folder(_rundir)
		_csFile = 	'%s/inputFiles/epig_input/chromStates/' \
					'%s/%s_chr%d_chromatin_states_From%dMbTo%dMb.txt'\
					%(glb_path,self.celltype,self.celltype,self.chrId,_gSta,_gEnd)
		_ctcfIndFile = 	'%s/inputFiles/epig_input/ctcfSites/' \
						'%s/%s_chr%d_ctcf_index_From%dMbTo%dMb.txt'\
						%(glb_path,self.celltype,self.celltype,self.chrId,_gSta,_gEnd)

		geninit = poly.GenInitConfig(_sepDist,self.lmpsdir,_ctcfFile)
		radius 	= geninit.confinement_size()

		# ----  lammps input file with custom random seed
		# ----  4928459 is simply a random seed number, could be changed to any other integer
		_inFile = _rundir + "in.chromosome"
		in_tmp = fi.input(self._lmpsTemplate)
		pf = open(_inFile, 'w')
		pf.write('variable        rseed equal   %d\n'%(4928459+self.runId) )			# 4928459 is simply a random seed number
																						# could be changed to any other integer
		for line in in_tmp:
			if line[0:9] == 'read_data':
				pf.write(
'''read_data       ../../../../inputFiles/lmps_input/%s/data.chromosome.chr%d\n'''%(self.celltype,self.chrId))
			elif line[0:10] == 'pair_style':
				if _sepDist == 25:
					ucs_name = 'ucs_chrom.txt'; uctcf_name = 'uctcf_chrom.txt';
				else:
					ucs_name = 'ucs_chrom_extension.txt'; uctcf_name = 'uctcf_chrom_extension.txt';
				pf.write(
'''pair_style      hybrid/overlay table linear 10000 tanhlr/cut/ideala 6.0 %d 15 %s/%s %s tanhlr/cut/ideal 6.0 %s/%s %s %d\n'''\
%(self._same_mol_flag_ideal,\
self._paramsFolder,ucs_name,_csFile,\
self._paramsFolder,uctcf_name,_ctcfIndFile,\
self.nearCtcfThreshold))
			elif line[0:16] == 'pair_coeff_softc':
				pf.write(
'''pair_coeff       * * table %s/soft_core_lj_4kT.table soft_core_lj 1.12\n'''%(self._paramsFolder))
			elif line[0:6] == 'region':
				pf.write('region             nucleus sphere 0.0 0.0 0.0 %.5f side in\n'%radius)
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
