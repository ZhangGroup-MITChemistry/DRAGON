from Ipt_module import *
from Params import *
Params()

import hoomd as hd
import hoomd.deprecated as hdp
import hoomd.md as hdmd

class GenInitConfig():

	_lmp_path = '%s/../src/md/lmps_input/'%glb_path
	_xml_path = '%s/initial_config/xml/'%_lmp_path
	_pdb_path = '%s/initial_config/pdb/'%_lmp_path
	_psf_path = '%s/initial_config/psf/'%_lmp_path
	_tcl_path = '%s/initial_config/tcl/'%_lmp_path
	_rdm_path = '%s/initial_config/lmps/rdm/'%_lmp_path
	_cds_path = '%s/initial_config/lmps/condense/'%_lmp_path

	def __init__(self,sepDist,lmpsdir):
		self.sepDist    = sepDist
		self.nbead      = int(sepDist*Mb/resolution)
		self.lmpsdir    = lmpsdir


	def genInitConfig(self):
		self.genxml()
		self.genpsfpdb()
		self.chrom_condense()
		self.move2center()


	def genxml(self):
		# ----  generate polymer configuration through random walk
		# ----  parameters
		phi_P = 0.05    # packing fraction
		n_poly = 1      # number of polymers
		T = 1.2         # temperature
		pType = ['A']*self.nbead
		hd.context.initialize("")
		polymer1 = dict(bond_len=1.2, type=pType,bond="linear", count=n_poly)
		
		# ----  find the length of the box
		N = len(polymer1['type']) * polymer1['count']

		# ----  generate the polymer system
		hdp.init.create_random_polymers(box=hd.data.boxdim(volume=pi*N/(6.0*phi_P)),
				polymers=[polymer1],separation=dict(A=0.35, B=0.35, C=0.35), seed=12)

		# ----  force field setup
		harmonic = hdmd.bond.harmonic()
		harmonic.bond_coeff.set('polymer', k=30.0, r0=1.0)
		nl = hdmd.nlist.cell()
		lj = hdmd.pair.lj(r_cut=3.0,nlist=nl)

		# ----  pure repulsive
		lj.pair_coeff.set('A', 'A', epsilon=1.0, sigma=1.0, alpha=1.0)
		all = hd.group.all()
		hdp.dump.xml(group=all,filename='%s/rdm_polymer_%d.xml'%(self._xml_path,self.sepDist),vis=True)
		hd.dump.dcd(filename='%s/rdm_polymer_%d.dcd'%(self._xml_path,self.sepDist),period=100,unwrap_full=True)
		
		# ----  integrate NVT for some time steps
		hdmd.integrate.mode_standard(dt=0.002)
		hdmd.integrate.nvt(group=all, kT=1.0, tau=0.5)
		hd.run(1000)


	def genpsfpdb(self):
		# ----  generate pdb and psf file
		tclScript='''
package require topotools

mol load hoomd   %s/rdm_polymer_%d.xml
animate read dcd %s/rdm_polymer_%d.dcd waitfor all

set nf [molinfo top get numframes]
puts $nf
animate write pdb %s/initial_config_%d.pdb beg [expr $nf-1]
animate write psf %s/initial_config_%d.psf
topo guessangles
topo writelammpsdata %s/data.chromosome_%d
'''%(self._xml_path,self.sepDist,\
self._xml_path,self.sepDist,\
self._pdb_path,self.sepDist,\
self._psf_path,self.sepDist,\
self._rdm_path,self.sepDist)

		fo = open('%s/genpsf_%d.tcl'%(self._tcl_path,self.sepDist), 'w')
		fo.write(tclScript)
		fo.close()

		cmd = 'vmd -dispdev text -eofexit < %s/genpsf_%d.tcl'%(self._tcl_path,self.sepDist)
		q = Popen(cmd, shell=True)
		q.communicate()


	def chrom_condense(self):
		# ----  condense the chromosome to fit into nucleus 
		inFile = open('%s/data.chromosome_%d'%(self._rdm_path,self.sepDist),'r')
		lines  = inFile.readlines()
		nl     = len(lines)

		fo = open('%s/data.chromosome'%self._cds_path,'w')
		for il in range(11):
			fo.write(lines[il])

		bd = self.nbead*1E2/resolution
		box = '''-%d   %d  xlo xhi
-%d   %d  ylo yhi
-%d   %d  zlo zhi
'''%(bd,bd,bd,bd,bd,bd)
		fo.write(box)

		for il in range(14, nl):
			fo.write(lines[il])
		fo.close()

		in_tmp = fi.input('%s/condense_lammps_template.in'%self._lmp_path)
		pf = open('%s/in.chromosome'%self._cds_path, 'w')
		for line in in_tmp:
			if line[0:7] == 'groupid':
				pf.write('group           chrom id 1:%d\n'%self.nbead)
			elif line[0:3] == 'run':
				pf.write('run           %d'%(self.nbead+5E3))
			else:
				pf.write(line)
		pf.close() 

		# ----	run MD
		cmd = 'cd %s; mpirun -np 1 %s/lmp_openmpi < in.chromosome;'%(self._cds_path,self.lmpsdir)
		q = Popen(cmd, shell=True)
		q.communicate()


	def move2center(self):
		# ----  move the chromosome to the center (0,0,0) for convenience
		# ----	the output would be the initial configuration for MD
		tclScript='''
package require topotools

mol load psf %s/initial_config_%d.psf pdb %s/initial_config_%d.pdb
animate read dcd %s/DUMP_FILE.dcd waitfor all
[atomselect top all] moveby [vecscale -1.0 [measure center [atomselect top all]]]

set sel [atomselect top all]
foreach sl [$sel get serial] {
	set sel [atomselect top "serial $sl"]
	$sel set name CA
	$sel set resname ALA
	$sel set resid $sl
	$sel delete 
}

set nf [molinfo top get numframes]
animate write pdb %s/initial_config_%d_centered.pdb beg [expr $nf-1]

#   ---
mol delete all
mol load psf %s/initial_config_%d.psf pdb %s/initial_config_%d_centered.pdb
topo guessangles
topo writelammpsdata %s/data.chromosome.init%d
'''%(self._psf_path,self.sepDist,self._pdb_path,self.sepDist,\
self._cds_path,\
self._pdb_path,self.sepDist,\
self._psf_path,self.sepDist,\
self._pdb_path,self.sepDist,\
self._lmp_path,self.sepDist)

		fo = open('%s/center_chromosome.tcl'%self._tcl_path, 'w')
		fo.write(tclScript)
		fo.close()

		cmd = 'vmd -dispdev text -eofexit < %s/center_chromosome.tcl'%self._tcl_path
		q = Popen(cmd, shell=True)
		q.communicate()


	def confinement_size(self):
		# ----  calculate the confinement size to have a consistent base pair density
		mb = self.nbead*resolution/Mb
		R = 50/(3.*(6.176/mb)**(1./3.))
		return R

