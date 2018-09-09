from Ipt_module import *

class Params():

	global Mb,resolution
	global glb_path,chr_region
	global comment,bond_coeffs,angle_coeffs
	global ncs,nctcf

#	---- default: 25Mb segment, at resolution of 5kb
	Mb=1E6
	resolution=5000

#	---- default: global path
	glb_path = os.getcwd()

#	---- chromosome segment region
	with open('%s/../src/chr_region.txt'%glb_path,'r') as f:
		chr_region = json.load(f)

#	---- LAMMPS input parameters
	comment = 'coor'
	bond_coeffs = [30.0, 1.5, 1.0, 1.0]
	angle_coeffs = [10.0, 30.0]

#	---- number of type of chromatin states and ctcfs
	ncs = 15
	nctcf = 4