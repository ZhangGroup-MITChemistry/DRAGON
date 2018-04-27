from Ipt_module import *

class Params():

	global Mb,resolution,nbead
	global glb_path,chr_region
	global comment,bond_coeffs,angle_coeffs
	global ncs,nctcf

#	---- default: 25Mb segment, at resolution of 5kb
	Mb=1E6
	resolution=5000
	nbead=int(25*Mb/resolution)

#	---- default: global path
	glb_path = os.getcwd()

#	---- chromosome segment region
	chr_region = np.loadtxt('%s/../src/chr_region.txt'%glb_path)

#	---- LAMMPS input parameters
	comment = 'coor'
	bond_coeffs = [30.0, 1.5, 1.0, 1.0]
	angle_coeffs = [10.0, 30.0]

#	---- number of type of chromatin states and ctcfs
	ncs = 15
	nctcf = 4