# generate the VMD scripts
from Ipt_module import *
from Params import *
Params()

import fileinput as fi

class VMD():

	def genPdb(self,celltype,chrId):

		gSta = chr_region[chrId-1,1]
		gEnd = gSta + 25
		cs_path = '%s/../../runMolecularDynamics/inputFiles/epig_input/chromStates/'\
			'/%s/%s_chr%d_chromatin_states_From%dMbTo%dMb.txt'\
			%(glb_path,celltype,celltype,chrId,gSta,gEnd)
		
		global pf_cs
		pf_cs = np.loadtxt(cs_path)
		global data_file_path
		data_file_path = '../../src/vmd/data_file'
		xyz = np.loadtxt('%s/atom_cor.txt'%data_file_path,usecols=(4,5,6))
		xyz = xyz.tolist()

		####
		pdb_path='./pdb/%s/'%celltype
		if not os.path.exists(pdb_path):
			os.makedirs(pdb_path)
		pf_1 = open('%s/chr%d_cs.pdb'%(pdb_path,chrId),'w')		# data file

		# Title
		pf_1.writelines('CRYST1   21.213   20.883   20.800  90.00  90.00  90.00 P 1           1\n')
		
		# Angle
		for i in range(1,nbead+1):
			pf_1.writelines('ATOM' + '%7d'%i + '%+4s'%('CA') +'  E%02d'%(pf_cs[i-1][1])+ ' X' + '%4d'%i + '%+12.3f'%(xyz[i-1][0]) + '%+8.3f'%(xyz[i-1][1]) + '%+8.3f'%(xyz[i-1][2]) +'\n')# + '  5.00  0.00      C\n')
		pf_1.writelines('END')
		pf_1.close()


	def genPsf(self,celltype,chrId):
		
		####
		psf_path='./psf/%s/'%celltype
		if not os.path.exists(psf_path):
			os.makedirs(psf_path)
		pf_1 = open('%s/chr%d_cs.psf'%(psf_path,chrId),'w')		# data file

		# Title
		pf_1.writelines('PSF\n\t1 !NTITLE\n REMARKS VMD-generated NAMD/X-Plor PSF structure file\n    5000 !NATOM\n')
		
		# Angle
		for i in range(1,nbead+1):
			pf_1.writelines('    ' + '%4d'%i + ' C    ' + '%4d'%i +' E%02d  '%(pf_cs[i-1][1])+ 'CA   CT     ' + '0.000000        1.0000           0\n')
		pf_1.write('\n')
		pf_1.close()

		f_out = open('%s/psf_template.txt'%data_file_path,'r')
		f_in = open('%s/chr%d_cs.psf'%(psf_path,chrId),'a+')
		for lines in f_out.readlines():
			f_in.writelines(lines)
		f_out.close()
		f_in.close()


	def genVMDScript(self,celltype,chrId):
		
		# ---- generate the VMD script ---- #
		# ---- visualize the individual chromosome ---- #

		vmd_tmp = fi.input('%s/VMDColorTemplate.vmd'%data_file_path)
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

