# generate the VMD scripts
from Ipt_module import *
from Params import *
Params()

class VMD():

	def modcs(self,pdbfilename,csmat):

		pdbDict = {0:'X',1:'A',2:'B',3:'C',4:'D',5:'E'}
		fi = open(pdbfilename)
		new = []
		count = 0
		for line in fi:
			if line[0:4] == 'ATOM':
				count += 1
				line = line.replace('ALA','E%02d'%csmat[count-1][1])
				idx = count//1E4
				line = line.replace('X',pdbDict[idx])
			new.append(line)
		return new

	def extractaddin(self,filename,nbead):

		# ----	this function returns addin lines from psf

		fi = open(filename)
		addin = []
		count = 0
		for line in fi:
			count += 1
			if count > nbead+7:
				addin.append(line)
		return addin

	def genPdb(self,celltype,chrId):

		# ----	this function could generate pdb file with chromatin states

		gSta = chr_region[chrId-1,1]
		gEnd = chr_region[chrId-1,2]
		global sepDist
		sepDist = gEnd-gSta
		global nbead
		nbead = int(sepDist*Mb/resolution)
		cs_path = '%s/../../runMolecularDynamics/inputFiles/epig_input/chromStates/'\
			'/%s/%s_chr%d_chromatin_states_From%dMbTo%dMb.txt'\
			%(glb_path,celltype,celltype,chrId,gSta,gEnd)
		
		global pf_cs
		pf_cs = np.loadtxt(cs_path)
		
		pdb_file_path = '../../src/md/lmps_input/initial_config/pdb/initial_config_%d_centered.pdb'%sepDist
		new = self.modcs(pdb_file_path,pf_cs)

		####
		pdb_path='./pdb/%s/'%celltype
		if not os.path.exists(pdb_path):
			os.makedirs(pdb_path)
		pf_1 = open('%s/chr%d_cs.pdb'%(pdb_path,chrId),'w')		# data file

		####
		for line in new:
			pf_1.writelines(line)
		pf_1.close()


	def genPsf(self,celltype,chrId):
		
		####
		global input_file_path
		input_file_path = '../../src/vmd/vmd_input'

		psf_path='./psf/%s/'%celltype
		if not os.path.exists(psf_path):
			os.makedirs(psf_path)
		pf_1 = open('%s/chr%d_cs.psf'%(psf_path,chrId),'w')		# data file

		# Title
		pf_1.writelines('PSF\n\n\t1 !NTITLE\n REMARKS VMD-generated NAMD/X-Plor PSF structure file\n\n   %d !NATOM\n'%nbead)
		
		# Angle
		for i in range(1,nbead+1):
			pf_1.writelines('   ' + '%5d'%i + ' C   ' + '%5d'%i +' E%02d  '%(pf_cs[i-1][1])+ 'CA   CT     ' + '0.000000        1.0000           0\n')
		pf_1.write('\n')
		pf_1.close()

		psf_file_path = '../../src/md/lmps_input/initial_config/psf/initial_config_%d.psf'%sepDist
		addin = self.extractaddin(psf_file_path,nbead)
		f_in = open('%s/chr%d_cs.psf'%(psf_path,chrId),'a+')
		for lines in addin:
			f_in.writelines(lines)
		f_in.close()


	def genVMDScript(self,celltype,chrId):
		
		# ---- generate the VMD script ---- #
		# ---- visualize the individual chromosome ---- #

		vmd_tmp = fi.input('%s/VMDColorTemplate.vmd'%input_file_path)
		vmd = open('./vmdScript/VMDColor_%s_chr%d.vmd'%(celltype,chrId),'w')
		for line in vmd_tmp:
			if line[0:19] == 'mol new chr1_cs.pdb':
				vmd.write('mol new ../pdb/%s/chr%d_cs.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'\
															%(celltype,chrId))
			elif line[0:23] == 'mol addfile chr1_cs.psf':
				vmd.write('mol addfile ../psf/%s/chr%d_cs.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'\
															%(celltype,chrId))
			elif line[0:26] == 'mol rename top chr1_cs.psf':
				vmd.write('mol rename top ../psf/%s/chr%d_cs.psf\n'%(celltype,chrId))
			else:
				vmd.write(line)
		vmd.close()
		print('''   > VMD script for for %s, chromosome %d is generated.'''\
					%(celltype,chrId))

