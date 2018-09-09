from Ipt_module import *

class Params():

	global Mb,resolution
	global cg_fac,num_kb
	global startb,endb
	global startfr
	global glb_path,chr_region
	
	# ---- default: 25Mb segment, at resolution of 5kb ---- #
	Mb=1E6
	resolution=5000

	# ---- default: global path
	glb_path = os.getcwd()

	#	---- chromosome segment region
	with open('%s/../../src/chr_region.txt'%glb_path,'r') as f:
		chr_region = json.load(f)
	
	# ---- default: visulize the map at 50kb ---- #
	cg_fac = 10
	num_kb = 5*cg_fac

	# ---- default: calculate the 20Mb (2Mb-22Mb) out of total 25Mb ---- #
	startb = 400
	endb = 4400

	# ---- default: start calculation from the 1st frame ---- #
	startfr = 1
	
	# ---- default: global path ---- #
	glb_path = os.getcwd()
	