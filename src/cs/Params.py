from Ipt_module import *

class Params():

	global Mb,resolution,nbead
	global glb_path
	global chr_region
	
	# ---- default: 25Mb segment, at resolution of 5kb ---- #
	Mb=1E6
	resolution=5000
	nbead=int(25*Mb/resolution)

	# ---- default: global path ---- #
	glb_path = os.getcwd()

	# ---- chromosome segment region ---- #
	chr_region = np.loadtxt('%s/../../../../src/chr_region.txt'%glb_path)
