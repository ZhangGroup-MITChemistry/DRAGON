from Ipt_module import *
from Params import *

import lammps_tools as lmp

class GenInitConfig():

    _lmp_path = '%s/../src/md/lmps_input/'%glb_path
    _pdb_path = '%s/initial_config/pdb/'%_lmp_path
    _psf_path = '%s/initial_config/psf/'%_lmp_path
    _tcl_path = '%s/initial_config/tcl/'%_lmp_path
    _rdm_path = '%s/initial_config/lmps/rdm/'%_lmp_path
    _cds_path = '%s/initial_config/lmps/condense/'%_lmp_path

    def __init__(self,sepDist,lmpsdir,ctcfFile):
        self.sepDist    = sepDist
        self.nbead      = int(sepDist*Mb/resolution)
        self.lmpsdir    = lmpsdir
        self.ctcfFile   = ctcfFile


    def genInitConfig(self):
        self.build_data()
        self.chrom_condense()
        self.genpsfpdb()
        self.move2center()


    def assignBond(self,hp):
    # ---- assign bond type
        hp.add_bond_type(bond_coeffs,comment)
        for atom in hp.atoms:
            ia = atom['i']
            if ia >= len(hp.atoms):
                break
            atoml = (ia, ia+1)
            hp.add_bond(atoml,comment=None,atom_names=None,i=None,bond_type=hp.bond_types[0])


    def assignAngle(self,hp):
    # ---- assign angle type
        hp.add_angle_type(angle_coeffs,comment)
        n_atom = len(hp.atoms)
        for atom in hp.atoms:        
            ia = atom['i']
            iap1 = ia+1
            if iap1 > n_atom:
                break
            iap2 = iap1+1
            if iap2 > n_atom:
                break
            atoml = (ia, iap1, iap2)
            hp.add_angle(atoml,comment=None,atom_names=None,i=None,angle_type=hp.angle_types[0])


    def assignAtom(self,hp):
    # ----  build xyz for a linear polymer
        xyz = np.zeros((self.nbead,3), dtype='float')
        xyz[:,1] = range(1,self.nbead+1)
        hp.add_atom_type(1.0, None, None)
        for ib in range(self.nbead):
            hp.add_atom(xyz[ib,0], xyz[ib,1], xyz[ib,2], 1, 0.0, comment=comment, i=ib+1, atom_type=hp.atom_types[0])


    def build_data(self):
    # ----  construct lammps data
        hp = lmp.Data()
        hp.header = 'Chromatin data file produced on: %s'%time.asctime()
        
        self.assignAtom(hp)
        self.assignBond(hp)
        self.assignAngle(hp)
        size=(-self.nbead,self.nbead);hp.box[0]=size;hp.box[1]=size;hp.box[2]=size;
        hp.write_to_file('%s/data.chromosome_%d'%(self._rdm_path,self.sepDist),ellipsoidFlag=0)


    def chrom_condense(self):
        # ----  condense the chromosome to fit into nucleus 
        inFile = open('%s/data.chromosome_%d'%(self._rdm_path,self.sepDist),'r')
        lines  = inFile.readlines()
        nl     = len(lines)

        fo = open('%s/data.chromosome'%self._cds_path,'w')
        for il in range(10):
            fo.write(lines[il])

        bd = self.nbead
        box = '''-%d   %d  xlo xhi
-%d   %d  ylo yhi
-%d   %d  zlo zhi
'''%(bd,bd,bd,bd,bd,bd)
        fo.write(box)

        for il in range(14, nl):
            fo.write(lines[il])
        inFile.close()
        fo.close()

        in_tmp = fi.input('%s/condense_lammps_template.in'%self._lmp_path)
        pf = open('%s/in.chromosome'%self._cds_path, 'w')
        for line in in_tmp:
            if line[0:7] == 'groupid':
                pf.write('group           chrom id 1:%d\n'%self.nbead)
            elif line[0:3] == 'run':
                pf.write('run           %d'%(self.nbead*5))
            else:
                pf.write(line)
        pf.close() 

        # ----  run MD
        cmd = 'cd %s; mpirun -np 1 %s/lmp_openmpi < in.chromosome;'%(self._cds_path,self.lmpsdir)
        q = Popen(cmd, shell=True)
        q.communicate()


    def genpsfpdb(self):
        # ----  generate pdb and psf file
        tclScript='''
package require topotools

topo readlammpsdata %s/data.chromosome_%d
animate read dcd %s/DUMP_FILE.dcd waitfor all

set nf [molinfo top get numframes]
puts $nf
animate write pdb %s/initial_config_%d.pdb beg [expr $nf-1]
animate write psf %s/initial_config_%d.psf
topo guessangles
'''%(self._rdm_path,self.sepDist,\
self._cds_path,\
self._pdb_path,self.sepDist,\
self._psf_path,self.sepDist)

        fo = open('%s/genpsf_%d.tcl'%(self._tcl_path,self.sepDist), 'w')
        fo.write(tclScript)
        fo.close()

        cmd = 'vmd -dispdev text -eofexit < %s/genpsf_%d.tcl'%(self._tcl_path,self.sepDist)
        q = Popen(cmd, shell=True)
        q.communicate()


    def move2center(self):
        # ----  move the chromosome to the center (0,0,0) for convenience
        # ----  the output would be the initial configuration for MD
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
self._psf_path,self.sepDist,self._pdb_path,self.sepDist,\
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

