from Ipt_module import *

class Params():

	global Mb,resolution
	global glb_path,chr_region

#	---- default: 25Mb segment, at resolution of 5kb
	Mb=1E6
	resolution=5000

#	---- default: global path
	glb_path = os.getcwd()

#	---- chromosome segment region
	with open('%s/../../src/chr_region.txt'%glb_path,'r') as f:
		chr_region = json.load(f)