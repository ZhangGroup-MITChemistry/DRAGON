from Ipt_module import *


class ProcessEpigData():

	global _glb_path
	_glb_path = os.path.abspath('.')

	def __init__(self,celltype,chrom_lst,bind_flxb,cap):

		self.celltype           = celltype
		self.chrom_lst          = chrom_lst
		self.bind_flxb 			= bind_flxb
		self.cap 				= cap


	def procChromState(self):
	# ----  
		os.chdir('%s/inputFiles/epig_input/chromStates/'%_glb_path)
		cmd = 'python genChromState.py -C %s -c '%self.celltype
		for chrId in self.chrom_lst:
			cmd += '%s '%str(chrId)
		os.system(cmd)


	def procCTCFSites(self):
	# ----  
		os.chdir('%s/inputFiles/epig_input/ctcfSites/'%_glb_path)
		cmd = 'python genCTCFbinding.py -C %s -b %d -a %d -c '\
						%(self.celltype,self.bind_flxb,self.cap)
		for chrId in self.chrom_lst:
			cmd += '%s '%str(chrId)
		os.system(cmd)
	