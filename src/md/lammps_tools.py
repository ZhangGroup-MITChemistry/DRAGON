# -*- coding: iso-8859-1 -*-
'''Python library to manage lammps files'''

import logging
import gzip
import itertools
import StringIO
import numpy as np
try:
  import pandas
except ImportError:
  pass

# import molecule

def load_file(filename='log.lammps', run=1):
  '''Helper function to load a lammps file'''

  f = open(filename, 'r')
  line = f.readline().strip()
  f.close()

  if line.startswith('# Time-averaged data for'):
    print 'Loading LAMMPS time-averaged data from', filename
    return Ave(filename).to_array()
  elif line.startswith('# Spatial-averaged data for'):
    print 'Loading LAMMPS spatial-average data from', filename
    return Ave(filename).to_array()
  elif line.startswith('LAMMPS'):
    print 'Loading LAMMPS log file from', filename, line[6:]
    return Log(filename, run).to_array()
  else:
    raise ValueError('Unknown format for file %s' % filename)

class Data(object):
  '''Class representing a Data file'''

  def __init__(self, filename=None):

    self.header = None

    self.atoms = []
    self.bonds = []
    self.angles = []
    self.dihedrals = []
    self.impropers = []
    self.extras = []

    self.atom_types = []
    self.bond_types = []
    self.angle_types = []
    self.dihedral_types = []
    self.improper_types = []

    self.mols = []

    self.box = [(None, None), (None, None), (None, None)]
    self.tilt = (0, 0, 0)

    if filename:
      self.read_from_file(filename)
      self.filename = filename
    else:
      self.filename = None

  def add_atom_type(self, mass, coeffs, comment=None, drude_type=None):
    '''Add an atom type to the end of the list'''

    i = len(self.atom_types) + 1
    atom_type = {'i': i, 'mass': mass, 'coeffs': coeffs, 'comment': comment, 'drude_type': drude_type}
    self.atom_types.append(atom_type)
    return i

  def add_bond_type(self, coeffs, comment=None, atom_names=None):
    '''Add a bond type to the end of the list'''

    if comment == 'auto':
      comment = '-'.join(atom_names)
    i = len(self.bond_types) + 1
    bond_type = {'i': i, 'atom_names': atom_names, 'coeffs': coeffs, 'comment': comment}
    self.bond_types.append(bond_type)
    return i

  def add_angle_type(self, coeffs, comment=None, atom_names=None):
    '''Add an angle type to the end of the list'''

    if comment == 'auto':
      comment = '-'.join(atom_names)
    i = len(self.angle_types) + 1
    angle_type = {'i': i, 'atom_names': atom_names, 'coeffs': coeffs, 'comment': comment}
    self.angle_types.append(angle_type)
    return i

  def add_dihedral_type(self, coeffs, comment=None, style=None, atom_names=None, is_improper=False):
    '''Add a dihedral type to the end of the list'''

    if comment == 'auto':
      comment = '-'.join(atom_names)
    i = len(self.dihedral_types) + 1
    dihedral_type = {'i': i, 'atom_names': atom_names, 'coeffs': coeffs, 'comment': comment, 'style': style, 'is_improper': is_improper}
    self.dihedral_types.append(dihedral_type)
    return i

  def add_improper_type(self, coeffs, comment=None, atom_names=None):
    '''Add a improper type to the end of the list'''

    if comment == 'auto':
      comment = '-'.join(atom_names)
    i = len(self.improper_types) + 1
    improper_type = {'i': i, 'atom_names': atom_names, 'coeffs': coeffs, 'comment': comment}
    self.improper_types.append(improper_type)
    return i

  def add_atom(self, x, y, z, mol_i, charge, mass=None, comment=None, i=None, atom_type=None):
    '''Add an atom in the list, by getting the type'''

    if atom_type is None:
      right_atom_type = None
      # Search for the right type
      for atom_type in self.atom_types:
        if atom_type['mass'] == mass and atom_type['comment'] == comment:
          right_atom_type = atom_type
          break

      if right_atom_type is None:
        raise ValueError('Atom type not found for atom %s' % comment)
    else:
      right_atom_type = atom_type

    if i is None:
      i = len(self.atoms) + 1
    atom = {'i': i, 'mol_i': mol_i, 'atom_type_i': right_atom_type['i'], 'atom_type': right_atom_type, \
            'charge': charge, 'x': x, 'y': y, 'z': z, 'comment': comment}
    try:
      self.atoms[i - 1] = atom
    except IndexError:
      # There are holes in the list, grow list (holes will be None and will be removed after)
      self.atoms.extend([None] * (i - len(self.atoms)))
      self.atoms[i - 1] = atom
    return i

  def add_bond(self, atom_is, comment=None, atom_names=None, i=None, bond_type=None):
    '''Add a bond in the list, by getting the type'''

    if bond_type is None:
      right_bond_type = None
      # Search for the right type
      for bond_type in self.bond_types:
        if bond_type['atom_names'] == atom_names:
          right_bond_type = bond_type
          is_reverse = False
          break
        if bond_type['atom_names'] == (atom_names[1], atom_names[0]):
          right_bond_type = bond_type
          is_reverse = True
          break
      if right_bond_type is None:
        raise ValueError('Bond type not found for bond %s' % comment)
    else:
      right_bond_type = bond_type

    if comment == 'auto':
      if is_reverse:
        comment = '-'.join(reversed(right_bond_type['atom_names']))
      else:
        comment = '-'.join(right_bond_type['atom_names'])

    atom1 = self.atoms[atom_is[0] - 1]
    atom2 = self.atoms[atom_is[1] - 1]

    if i is None:
      i = len(self.bonds) + 1
    bond = {'i': i, 'bond_type_i': right_bond_type['i'], 'bond_type': right_bond_type, \
            'atom1': atom1, 'atom2': atom2, 'comment': comment}
    try:
      self.bonds[i - 1] = bond
    except IndexError:
      # There are holes in the list, grow list (holes will be None and will be removed after)
      self.bonds.extend([None] * (i - len(self.bonds)))
      self.bonds[i - 1] = bond
    return i

  def add_angle(self, atom_is, comment=None, atom_names=None, i=None, angle_type=None):
    '''Add an angle in the list, by getting the type'''

    if angle_type is None:
      right_angle_type = None
      # Search for the right type
      for angle_type in self.angle_types:
        if angle_type['atom_names'] == atom_names:
          right_angle_type = angle_type
          is_reverse = False
          break
        if angle_type['atom_names'] == (atom_names[2], atom_names[1], atom_names[0]):
          right_angle_type = angle_type
          is_reverse = True
          break

      if right_angle_type is None:
        raise ValueError('Angle type not found for angle %s' % comment)
    else:
      right_angle_type = angle_type

    atom1 = self.atoms[atom_is[0] - 1]
    atom2 = self.atoms[atom_is[1] - 1]
    atom3 = self.atoms[atom_is[2] - 1]

    if comment == 'auto':
      if is_reverse:
        comment = '-'.join(reversed(right_angle_type['atom_names']))
      else:
        comment = '-'.join(right_angle_type['atom_names'])

    if i is None:
      i = len(self.angles) + 1
    angle = {'i': i, 'angle_type_i': right_angle_type['i'], 'angle_type': right_angle_type, \
             'atom1': atom1, 'atom2': atom2, 'atom3': atom3, 'comment': comment}
    try:
      self.angles[i - 1] = angle
    except IndexError:
      # There are holes in the list, grow list (holes will be None and will be removed after)
      self.angles.extend([None] * (i - len(self.angles)))
      self.angles[i - 1] = angle
    return i

  def add_dihedral(self, atom_is, comment=None, atom_names=None, i=None, dihedral_type=None, is_improper=False):
    '''Add an dihedral in the list, by getting the type'''

    if dihedral_type is None:
      right_dihedral_type = None
      # Search for the right type
      for dihedral_type in self.dihedral_types:
        if dihedral_type['is_improper'] == is_improper:
          if dihedral_type['atom_names'] == atom_names:
            right_dihedral_type = dihedral_type
            is_reverse = False
            break
          if dihedral_type['atom_names'] == (atom_names[3], atom_names[2], atom_names[1], atom_names[0]):
            right_dihedral_type = dihedral_type
            is_reverse = True
            break

      if right_dihedral_type is None:
        raise ValueError('Dihedral type not found for dihedral %s' % comment)
    else:
      right_dihedral_type = dihedral_type

    atom1 = self.atoms[atom_is[0] - 1]
    atom2 = self.atoms[atom_is[1] - 1]
    atom3 = self.atoms[atom_is[2] - 1]
    atom4 = self.atoms[atom_is[3] - 1]

    if comment == 'auto':
      if is_reverse:
        comment = '-'.join(reversed(right_dihedral_type['atom_names']))
      else:
        comment = '-'.join(right_dihedral_type['atom_names'])

    if i is None:
      i = len(self.dihedrals) + 1
    dihedral = {'i': i, 'dihedral_type_i': right_dihedral_type['i'], 'dihedral_type': right_dihedral_type, \
                'atom1': atom1, 'atom2': atom2, 'atom3': atom3, 'atom4': atom4, 'comment': comment}
    try:
      self.dihedrals[i - 1] = dihedral
    except IndexError:
      # There are holes in the list, grow list (holes will be None and will be removed after)
      self.dihedrals.extend([None] * (i - len(self.dihedrals)))
      self.dihedrals[i - 1] = dihedral
    return i

  def add_improper(self, atom_is, comment=None, atom_names=None, i=None, improper_type=None):
    '''Add an improper in the list, by getting the type'''

    if improper_type is None:
      right_improper_type = None
      # Search for the right type
      for improper_type in self.improper_types:
        if improper_type['atom_names'] == atom_names or improper_type['atom_names'] == (atom_names[2], atom_names[1], atom_names[0]):
          right_improper_type = improper_type
          break

      if right_improper_type is None:
        raise ValueError('Improper type not found for improper %s' % comment)
    else:
      right_improper_type = improper_type

    atom1 = self.atoms[atom_is[0] - 1]
    atom2 = self.atoms[atom_is[1] - 1]
    atom3 = self.atoms[atom_is[2] - 1]
    atom4 = self.atoms[atom_is[3] - 1]

    if comment == 'auto':
      comment = '-'.join(right_improper_type['atom_names'])

    if i is None:
      i = len(self.impropers) + 1
    improper = {'i': i, 'improper_type_i': right_improper_type['i'], 'improper_type': right_improper_type, \
                'atom1': atom1, 'atom2': atom2, 'atom3': atom3, 'atom4': atom4, 'comment': comment}
    try:
      self.impropers[i - 1] = improper
    except IndexError:
      # There are holes in the list, grow list (holes will be None and will be removed after)
      self.impropers.extend([None] * (i - len(self.impropers)))
      self.impropers[i - 1] = improper
    return i

  def add_extra(self, fmt, values, comment=None, i=None):
    '''Add some infor for the extra part'''

    tokens = []
    for (f, value) in zip(fmt, values):
      tokens.append(f % value)

    if i is None:
      i = len(self.extras) + 1
    extra = {'i': i, 'line': ' '.join(tokens), 'comment': comment}
    try:
      self.extras[i - 1] = extra
    except IndexError:
      # There are holes in the list, grow list (holes will be None and will be removed after)
      self.extras.extend([None] * (i - len(self.extras)))
      self.extras[i - 1] = extra
    return i

  def read_from_file(self, filename='lammps.data', ellipsoidFlag=0):

    nb_atoms = 0
    nb_bonds = 0
    nb_angles = 0
    nb_dihedrals = 0
    nb_impropers = 0
    nb_atom_types = 0
    nb_bond_types = 0
    nb_angle_types = 0
    nb_dihedral_types = 0
    nb_improper_types = 0

    nb_mols = 0

    f = open(filename, 'r')
    self.header = f.readline().strip()

    section = 'Header'

    for line in f:
      line = line.strip()

      if line == '':
        continue

      keyword = line.split('#')[0].strip()

      # Section change
      if keyword == 'Masses' or keyword == 'Pair Coeffs' or keyword == 'Bond Coeffs' or keyword == 'Angle Coeffs' or keyword == 'Dihedral Coeffs' or keyword == 'Improper Coeffs' or \
         keyword == 'Atoms' or keyword == 'Bonds' or keyword == 'Angles' or keyword == 'Dihedrals' or keyword == 'Impropers' or keyword == 'Velocities' or \
         keyword == 'Drude Types' or keyword == 'Extras' or keyword == 'EXTRA' or \
         keyword == 'Ellipsoids':
        section = keyword
        continue

      # Header
      if section == 'Header':
        if 'atoms' in line:
          nb_atoms = int(line.split()[0])
          self.atoms = [None] * nb_atoms
          continue
        if 'bonds' in line:
          nb_bonds = int(line.split()[0])
          self.bonds = [None] * nb_bonds
          continue
        if 'angles' in line:
          nb_angles = int(line.split()[0])
          self.angles = [None] * nb_angles
          continue
        if 'dihedrals' in line:
          nb_dihedrals = int(line.split()[0])
          self.dihedrals = [None] * nb_dihedrals
          continue
        if 'impropers' in line:
          nb_impropers = int(line.split()[0])
          self.impropers = [None] * nb_impropers
          continue

        if 'atom types' in line:
          nb_atom_types = int(line.split()[0])
          continue
        if 'bond types' in line:
          nb_bond_types = int(line.split()[0])
          continue
        if 'angle types' in line:
          nb_angle_types = int(line.split()[0])
          continue
        if 'dihedral types' in line:
          nb_dihedral_types = int(line.split()[0])
          continue
        if 'improper types' in line:
          nb_improper_types = int(line.split()[0])
          continue

        if 'xlo xhi' in line:
          self.box[0] = (float(line.split()[0]), float(line.split()[1]))
          continue
        if 'ylo yhi' in line:
          self.box[1] = (float(line.split()[0]), float(line.split()[1]))
          continue
        if 'zlo zhi' in line:
          self.box[2] = (float(line.split()[0]), float(line.split()[1]))
          continue
        if 'xy xz yz' in line:
          self.tilt = (float(line.split()[0]), float(line.split()[1]), float(line.split()[2]))
          continue

      # Masses
      if section == 'Masses':
        if line.split()[0] == '#':
          continue
        i = int(line.split()[0])
        mass = float(line.split()[1])

        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None
        # Coeffs are not set yet
        self.add_atom_type(mass, None, comment)

      # Pair Coeffs
      if section == 'Pair Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break

        atom_type = self.atom_types[i - 1]
        atom_type['coeffs'] = coeffs

      # Drude Types
      if section == 'Drude Types':
        i = int(line.split()[0])
        # Get all float coeffs
        try:
          drude_type = int(line.split()[1])
        except ValueError:
          break

        atom_type = self.atom_types[i - 1]
        atom_type['drude type'] = drude_type

      # Bond Coeffs
      if section == 'Bond Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None
        self.add_bond_type(coeffs, comment)

      # Angle Coeffs
      if section == 'Angle Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None
        self.add_angle_type(coeffs, comment)

      # Dihedral Coeffs
      if section == 'Dihedral Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None
        self.add_dihedral_type(coeffs, comment)

      # Improper Coeffs
      if section == 'Improper Coeffs':
        i = int(line.split()[0])
        # Get all float coeffs
        coeffs = []
        for c in line.split()[1:]:
          try:
            coeff = float(c)
            coeffs.append(coeff)
          except ValueError:
            break
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None
        self.add_improper_type(coeffs, comment)

      # Atoms
      #TODO assume atom_style is full
      if section == 'Atoms':
        if line.split()[0] == '#':
          continue
        if ellipsoidFlag == 0:
          i = int(line.split()[0])
          mol_i = int(line.split()[1])
          atom_type_i = int(line.split()[2])
          charge = float(line.split()[3])
          x = float(line.split()[4])
          y = float(line.split()[5])
          z = float(line.split()[6])
          if '#' in line:
            comment = line.split('#')[1].strip()
          else:
            comment = None

          self.add_atom(x, y, z, mol_i, charge, comment=comment, i=i, atom_type=self.atom_types[atom_type_i - 1])
          if mol_i > nb_mols:
            nb_mols = mol_i
        else:
          i = int(line.split()[0])
          atom_type_i = int(line.split()[1])
          x = float(line.split()[2])
          y = float(line.split()[3])
          z = float(line.split()[4])

          mol_i = int(line.split()[5])
          charge = float(line.split()[6])

          if '#' in line:
            comment = line.split('#')[1].strip()
          else:
            comment = None

          self.add_atom(x, y, z, mol_i, charge, comment=comment, i=i, atom_type=self.atom_types[atom_type_i - 1])
          if mol_i > nb_mols:
            nb_mols = mol_i

      # Ellipsoids
      #TODO Not using any of the information here
      if section == 'Ellipsoids':
        continue

      # Bonds
      if section == 'Bonds':
        i = int(line.split()[0])
        bond_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None
        
        while len(self.bond_types) < bond_type_i:
          coeffs = []
          comment = None
          self.add_bond_type(coeffs, comment)

        self.add_bond((atom1_i, atom2_i), comment=comment, i=i, bond_type=self.bond_types[bond_type_i - 1])

      # Angles
      if section == 'Angles':
        i = int(line.split()[0])
        angle_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])
        atom3_i = int(line.split()[4])
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None

        while len(self.angle_types) < angle_type_i:
          coeffs = []
          comment = None
          self.add_angle_type(coeffs, comment)

        self.add_angle((atom1_i, atom2_i, atom3_i), comment=comment, i=i, angle_type=self.angle_types[angle_type_i - 1])

      # Dihedrals
      if section == 'Dihedrals':
        i = int(line.split()[0])
        dihedral_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])
        atom3_i = int(line.split()[4])
        atom4_i = int(line.split()[5])
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None

        self.add_dihedral((atom1_i, atom2_i, atom3_i, atom4_i), comment=comment, i=i, dihedral_type=self.dihedral_types[dihedral_type_i - 1])

      # Impropers
      if section == 'Impropers':
        i = int(line.split()[0])
        improper_type_i = int(line.split()[1])
        atom1_i = int(line.split()[2])
        atom2_i = int(line.split()[3])
        atom3_i = int(line.split()[4])
        atom4_i = int(line.split()[5])
        if '#' in line:
          comment = line.split('#')[1].strip()
        else:
          comment = None

        self.add_improper((atom1_i, atom2_i, atom3_i, atom4_i), comment=comment, i=i, improper_type=self.improper_types[improper_type_i - 1])

      # Velocities
      if section == 'Velocities':
        i = int(line.split()[0])
        vx = float(line.split()[1])
        vy = float(line.split()[2])
        vz = float(line.split()[3])
        self.atoms[i-1]['velocities'] = (vx, vy, vz)

      # Extra
      if section == 'EXTRA' or section == 'Extras':
        i = int(line.split()[0])
        # Get extra fields
        extra = line.split(None, 1)[1]
        if '#' in line:
          extra = extra.split('#')[0].strip()
        extras = [float(e) for e in extra.split()]
        self.atoms[i-1]['extras'] = tuple(extras)

    f.close()

    # Remove holes (None)
    if None in self.atoms:
      tmp_atoms = self.atoms
      self.atoms = [atom for atom in tmp_atoms if atom is not None]
    if None in self.bonds:
      tmp_bonds = self.bonds
      self.bonds = [bond for bond in tmp_bonds if bond is not None]
    if None in self.angles:
      tmp_angles = self.angles
      self.angles = [angle for angle in tmp_angles if angle is not None]
    if None in self.dihedrals:
      tmp_dihedrals = self.dihedrals
      self.dihedrals = [dihedral for dihedral in tmp_dihedrals if dihedral is not None]
    if None in self.impropers:
      tmp_impropers = self.impropers
      self.impropers = [improper for improper in tmp_impropers if improper is not None]
    if None in self.extras:
      tmp_extras = self.extras
      self.extras = [extra for extra in tmp_extras if extra is not None]

    # Validate nbs
    if len(self.atoms) != nb_atoms:
      raise ValueError('Nb of atoms is not coherent')
    if len(self.bonds) != nb_bonds:
      raise ValueError('Nb of bonds is not coherent')
    if len(self.angles) != nb_angles:
      raise ValueError('Nb of angles is not coherent')
    if len(self.dihedrals) != nb_dihedrals:
      raise ValueError('Nb of dihedrals is not coherent')
    if len(self.impropers) != nb_impropers:
      raise ValueError('Nb of impropers is not coherent')

    if len(self.atom_types) != nb_atom_types:
      raise ValueError('Nb of atom types is not coherent')
    if len(self.bond_types) != nb_bond_types:
      raise ValueError('Nb of bond types is not coherent')
    if len(self.angle_types) != nb_angle_types:
      raise ValueError('Nb of angle types is not coherent')
    if len(self.dihedral_types) != nb_dihedral_types:
      raise ValueError('Nb of dihedral types is not coherent')
    if len(self.improper_types) != nb_improper_types:
      raise ValueError('Nb of improper types is not coherent')

    # Reconstruct mols
    self.mols = []
    for i in xrange(nb_mols):
      mol = {'i': i + 1}
      self.mols.append(mol)

    for atom in self.atoms:
      mol = self.mols[atom['mol_i'] - 1]
      if 'atoms' not in mol:
        mol['atoms'] = []
      mol['atoms'].append(atom)
      atom['mol'] = mol

  def write_to_file(self, filename='lammps.data', ellipsoidFlag=0, coeffFlag=0):
    '''Write the data to a file'''

    # For stdout writing
    if isinstance(filename, file):
      f = filename
    else:
      f = open(filename, 'w')

    # Header
    f.write('%s\n\n' % self.header)
    f.write('%d atoms\n' % len(self.atoms))
    if ellipsoidFlag ==1:
      f.write('%d ellipsoids\n' % len(self.atoms))
    if len(self.bonds) > 0:
      f.write('%d bonds\n' % len(self.bonds))
    if len(self.angles) > 0:
      f.write('%d angles\n' % len(self.angles))
    if len(self.dihedrals) > 0:
      f.write('%d dihedrals\n' % len(self.dihedrals))
    if len(self.impropers) > 0:
      f.write('%d impropers\n' % len(self.impropers))
    f.write('\n')

    f.write('%d atom types\n' % len(self.atom_types))
    if len(self.bond_types) > 0:
      f.write('%d bond types\n' % len(self.bond_types))
    if len(self.angle_types) > 0:
      f.write('%d angle types\n' % len(self.angle_types))
    if len(self.dihedral_types) > 0:
      f.write('%d dihedral types\n' % len(self.dihedral_types))
    if len(self.improper_types) > 0:
      f.write('%d improper types\n' % len(self.improper_types))
    f.write('\n')

    f.write('%f %f xlo xhi\n' % self.box[0])
    f.write('%f %f ylo yhi\n' % self.box[1])
    f.write('%f %f zlo zhi\n' % self.box[2])
    if self.tilt != (0, 0, 0):
      f.write('%f %f %f xy xz yz\n' % self.tilt)
    f.write('\n')

    # Masses
    f.write('Masses\n\n')
    i = 0
    for atom_type in self.atom_types:
      i += 1
      if 'comment' in atom_type and atom_type['comment'] is not None:
        f.write('%d %d # %s\n' % (i, atom_type['mass'], atom_type['comment']))
      else:
        f.write('%d %d\n' % (i, atom_type['mass']))
    f.write('\n')

    # Pair Coeffs
    # Check if the pair_coeffs have to be written here
    if None not in [atom_type['coeffs'] for atom_type in self.atom_types]:
      f.write('Pair Coeffs\n\n')
      i = 0
      for atom_type in self.atom_types:
        i += 1
        coeffs_str = ' '.join(['%9.4f' % coeff for coeff in atom_type['coeffs']])
        if 'comment' in atom_type and atom_type['comment'] is not None:
          f.write('%4d %s # %s\n' % (i, coeffs_str, atom_type['comment']))
        else:
          f.write('%4d %s\n' % (i, coeffs_str))
      f.write('\n')

    # Drude Types
    if None not in [atom_type['drude_type'] for atom_type in self.atom_types]:
      f.write('Drude Types\n\n')
      i = 0
      for atom_type in self.atom_types:
        i += 1
        if 'comment' in atom_type and atom_type['comment'] is not None:
          f.write('%4d %4d # %s\n' % (i, atom_type['drude_type'], atom_type['comment']))
        else:
          f.write('%4d %4d\n' % (i, atom_type['drude_type']))
      f.write('\n')

    # Bond Coeffs
    if (self.bond_types != []) and (coeffFlag == 1):
      f.write('Bond Coeffs\n\n')
      i = 0
      for bond_type in self.bond_types:
        i += 1
        coeffs_str = ' '.join(['%9.4f' % coeff for coeff in bond_type['coeffs']])
        if 'comment' in bond_type and bond_type['comment'] is not None:
          f.write('%4d %s # %s\n' % (i, coeffs_str, bond_type['comment']))
        else:
          f.write('%4d %s\n' % (i, coeffs_str))
      f.write('\n')

    # Angle Coeffs
    if (self.angle_types != []) and (coeffFlag == 1):
      f.write('Angle Coeffs\n\n')
      i = 0
      for angle_type in self.angle_types:
        i += 1
        coeffs_str = ' '.join(['%9.4f' % coeff for coeff in angle_type['coeffs']])
        if 'comment' in angle_type and angle_type['comment'] is not None:
          f.write('%4d %s # %s\n' % (i, coeffs_str, angle_type['comment']))
        else:
          f.write('%4d %s\n' % (i, coeffs_str))
      f.write('\n')

    # Dihedral Coeffs
    if (self.dihedral_types != []) and (coeffFlag == 1):
      # If there's no style defined, there's only one
      try:
        nb_style = len(set([d['style'] for d in self.dihedral_types]))
      except KeyError:
        nb_style = 1
      f.write('Dihedral Coeffs\n\n')
      i = 0
      for dihedral_type in self.dihedral_types:
        i += 1
        coeffs_str = ' '.join(['%9.4f' % coeff for coeff in dihedral_type['coeffs']])
        # Don't forget the style if multiple
        if nb_style == 1:
          if 'comment' in dihedral_type and dihedral_type['comment'] is not None:
            f.write('%4d %s # %s\n' % (i, coeffs_str, dihedral_type['comment']))
          else:
            f.write('%4d %s\n' % (i, coeffs_str))
        else:
          if 'comment' in dihedral_type and dihedral_type['comment'] is not None:
            f.write('%4d %7s %s # %s\n' % (i, dihedral_type['style'], coeffs_str, dihedral_type['comment']))
          else:
            f.write('%4d %7s %s\n' % (i, dihedral_type['style'], coeffs_str))
      f.write('\n')

    # Improper Coeffs
    if (self.improper_types != []) and (coeffFlag == 1):
      f.write('Improper Coeffs\n\n')
      i = 0
      for improper_type in self.improper_types:
        i += 1
        coeffs_str = ' '.join(['%9.4f' % coeff for coeff in improper_type['coeffs']])
        if 'comment' in improper_type and improper_type['comment'] is not None:
          f.write('%4d %s # %s\n' % (i, coeffs_str, improper_type['comment']))
        else:
          f.write('%4d %s\n' % (i, coeffs_str))
      f.write('\n')

    # Atoms
    if ellipsoidFlag == 0:
        f.write('Atoms\n\n')
        for atom in self.atoms:
            #TODO assume atom_style is full
            if 'comment' in atom and atom['comment'] is not None:
                f.write('%d %d %d %.4f %13.6e %13.6e %13.6e # %s\n' % (atom['i'], atom['mol_i'], atom['atom_type_i'], atom['charge'], atom['x'], atom['y'], atom['z'], atom['comment']))
                # f.write('%d %d %f %f %f %d 0 1.90986 # %s\n' % (atom['i'], atom['atom_type_i'], atom['x'], atom['y'], atom['z'], atom['mol_i'], atom['comment']))
            else:
                #f.write('%d %d %d %.4f %13.6e %13.6e %13.6e\n' % (atom['i'], atom['mol_i'], atom['atom_type_i'], atom['charge'], atom['x'], atom['y'], atom['z']))
                f.write('%d %d %d %f %f %f 1 0 1.90986\n' % (atom['i'], atom['mol_i'], atom['atom_type_i'], atom['x'], atom['y'], atom['z']))
    else:
        f.write('# Atom-ID, type, position, molecule-ID, charge, ellipsoid flag, density\n')
        f.write('Atoms\n\n')
        # write the hybrid format of full and ellipsoid
        for atom in self.atoms:
            f.write('%d %d %f %f %f %d 1 1.90986\n' % (atom['i'], atom['atom_type_i'], atom['x'], atom['y'], atom['z'], atom['mol_i']))#, atom['charge']))

    f.write('\n')

    if ellipsoidFlag == 1:
        f.write('# Atom-ID, shape, quaternion\n')
        f.write('Ellipsoids\n')
        f.write('\n')
        for atom in self.atoms:
            f.write('%d 1 1 1 1 0 0 0\n'%(atom['i']))
        f.write('\n')

    # Bonds
    if self.bonds != []:
      f.write('Bonds\n\n')
      for bond in self.bonds:
        if 'comment' in bond and bond['comment'] is not None:
          f.write('%d %d %d %d # %s\n' % (bond['i'], bond['bond_type']['i'], bond['atom1']['i'], bond['atom2']['i'], bond['comment']))
        else:
          f.write('%d %d %d %d\n' % (bond['i'], bond['bond_type']['i'], bond['atom1']['i'], bond['atom2']['i']))
      f.write('\n')

    # Angles
    if self.angles != []:
      f.write('Angles\n\n')
      for angle in self.angles:
        if 'comment' in angle and angle['comment'] is not None:
          f.write('%d %d %d %d %d # %s\n' % (angle['i'], angle['angle_type']['i'], angle['atom1']['i'], angle['atom2']['i'], angle['atom3']['i'], angle['comment']))
        else:
          f.write('%d %d %d %d %d\n' % (angle['i'], angle['angle_type']['i'], angle['atom1']['i'], angle['atom2']['i'], angle['atom3']['i']))
      f.write('\n')

    # Dihedrals
    if self.dihedrals != []:
      f.write('Dihedrals\n\n')
      for dihedral in self.dihedrals:
        if 'comment' in dihedral and dihedral['comment'] is not None:
          f.write('%7d %7d %7d %7d %7d %7d # %s\n' % (dihedral['i'], dihedral['dihedral_type']['i'], \
                  dihedral['atom1']['i'], dihedral['atom2']['i'], dihedral['atom3']['i'], dihedral['atom4']['i'], dihedral['comment']))
        else:
          f.write('%7d %7d %7d %7d %7d %7d\n' % (dihedral['i'], dihedral['dihedral_type']['i'], dihedral['atom1']['i'], dihedral['atom2']['i'], dihedral['atom3']['i'], dihedral['atom4']['i']))
      f.write('\n')

    # Impropers
    if self.impropers != []:
      f.write('Impropers\n\n')
      for improper in self.impropers:
        if 'comment' in improper and improper['comment'] is not None:
          f.write('%7d %7d %7d %7d %7d %7d # %s\n' % (improper['i'], improper['improper_type']['i'], \
                  improper['atom1']['i'], improper['atom2']['i'], improper['atom3']['i'], improper['atom4']['i'], improper['comment']))
        else:
          f.write('%7d %7d %7d %7d %7d %7d\n' % (improper['i'], improper['improper_type']['i'], improper['atom1']['i'], improper['atom2']['i'], improper['atom3']['i'], improper['atom4']['i']))

    # Extra
    if self.extras != []:
      f.write('Extras\n\n')
      for extra in self.extras:
        if 'comment' in extra and extra['comment'] is not None:
          f.write('%7d %s # %s\n' % (extra['i'], extra['line'], extra['comment']))
        else:
          f.write('%7d %s\n' % (extra['i'], extra['line']))

    # Velocities

    # TODO ignore Velocities ?

class Log(object):
  '''Class representing a Log file reader (for thermo infos)'''

  def __init__(self, filename='log.lammps', run=1):

    self.filename = filename

    self.run = run

    # Default is read only one run
    self.nb_runs = 1

    self.last_line = False

    self.steps = 0

    self.is_multi = False

    # Temp values
    self.next_step = None
    self.next_cpu = None

    self.fields = []

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = gzip.open(filename, 'r')
    except IOError:
      self.f = open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = open(filename, 'r')

    # Seek to the start of the thermo stuff (and skip)
    for _ in xrange(self.run):
      self.seek_next_run()
    if self.run == 0:
      self.seek_next_run()

  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    stats = {}

    if self.run == 0:
      stats['run'] = self.nb_runs
    else:
      stats['run'] = self.run

    if not self.is_multi:

      line = self.f.readline()
      if line == '':
        raise StopIteration
      if line.startswith('Loop'):
        if self.run == 0:
          self.seek_next_run(True)
          return self.next()
        else:
          raise StopIteration
      values = line.split()

      # Skip SHAKE lines
      while 'SHAKE' in line or len(values) != len(self.fields):
        line = self.f.readline()
        if line == '':
          raise StopIteration
        if line.startswith('Loop'):
          if self.run == 0:
            self.seek_next_run(True)
            return self.next()
          else:
            raise StopIteration
        values = line.split()

      for (field, value) in zip(self.fields, values):
        # Special case for int values
        if field in ('step', 'elapsed'):
          stats[field] = int(value)
        else:
          stats[field] = float(value)
    else:
      if self.last_line:
        raise StopIteration

      # Multi lines version
      stats['step'] = self.next_step
      stats['cpu'] = self.next_cpu

      line = self.f.readline()
      if line == '':
        raise StopIteration
      values = line.split()
      # Get all the data in format: Field = Value
      while len(values) > 2 and values[1] == '=':
        # For each line, get the infos (line format: Field = Value)
        field = None
        for value in values:
          if value == '=':
            continue
          if field is None:
            field = value
          else:
            stats[field.lower()] = float(value)
            field = None
        line = self.f.readline()
        values = line.split()

      # Remove volume for NVT
      if 'volume' in self.fields and 'volume' not in stats:
        self.fields.remove('volume')

      # Try to get the next record
      while not line.startswith('-----') and line != '' and not line.startswith('Memory usage per processor'):
        line = self.f.readline()

      # Reach next log
      if line.startswith('Memory usage per processor'):
        if self.run == 0:
          self.seek_next_run(True, line=line)
        else:
          # It's over, but at next next
          self.last_line = True

      if line == '':
        self.last_line = True

      # We reach the first line of the next record
      if line.startswith('-----'):
        values = line.split()
        try:
          self.next_step = int(values[2])
          self.next_cpu = float(values[6])
        except IndexError:
          self.last_line = True

    return stats

  def to_array(self):
    '''Convert the log to a numpy array'''

    try:
      nb_fields = len(self.fields)
      types = zip(self.fields, [float] * nb_fields)
      return np.array([tuple([line[field] for field in self.fields]) for line in self], dtype=types)
    except ValueError:
      # It's NVT in fact, retry, volume has been already removed !
      nb_fields = len(self.fields)
      types = zip(self.fields, [float] * nb_fields)
      return np.array([tuple([line[field] for field in self.fields]) for line in self], dtype=types)

  def seek_next_run(self, is_multi_run=False, line=None):
    '''Seek to next run'''

    if not line:
      line = self.f.readline()
    while not line.startswith('Memory usage per processor'):
      line = self.f.readline()
      # Get and store steps (if int)
      if line.startswith('run'):
        try:
          self.steps = int(line.split()[1])
        except ValueError:
          pass
      # If EOF
      if line == '':
        if is_multi_run:
          raise StopIteration
        elif self.run != 1 and self.run != 0:
          raise ValueError('File %s does not seem to be a valid LAMMPS log file (with thermo infos) or does not contain info for run %d' % (self.filename, self.run))
        else:
          raise ValueError('File %s does not seem to be a valid LAMMPS log file (with thermo infos)' % self.filename)

    if is_multi_run:
      self.nb_runs += 1

    # Read the header to store in fields
    line = self.f.readline()
    if not line.startswith('-----'):
      tmp_fields = [field.lower() for field in line.split()]
      self.fields = []
      mult = {}
      for field in tmp_fields:
        # Manage multiple fields with same name
        if tmp_fields.count(field) > 1:
          if field not in mult:
            mult[field] = 0
          self.fields.append('%s_%d' % (field, mult[field] + 1))
          mult[field] += 1
        else:
          self.fields.append(field)
    else:
      self.next_step = int(line.split()[2])
      self.next_cpu = float(line.split()[6])
      self.is_multi = True
      # Suppose it's NPT (volume is here), we'll remove it after
      self.fields = ['step', 'cpu', 'toteng', 'kineng', 'temp', 'poteng', 'e_bond', 'e_angle', 'e_dihed', 'e_impro', 'e_vdwl', 'e_coul', 'e_long', 'press', 'volume']


class Dump(object):
  '''Class representing a Dump file reader'''

  def __init__(self, filename='lammps.dump', skip=0, sort=True, numpy=False, data=None, raw=True):

    self.filename = filename

    self.sort = sort
    self.numpy = numpy
    self.data = data
    self.raw = raw
    self.i = skip + 1
    self.atom_types = {}

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Test the read
      line = self.f.readline().strip()
      if line != 'ITEM: TIMESTEP':
        raise IOError('Not a valid LAMMPS dump file')
      self.f.readline()
      self.f.readline()
      # And get the nb of atoms
      nbatoms = int(self.f.readline())
      self.f.close()
      # If OK, reopen it
      self.f = gzip.open(filename, 'r')
    except IOError:
      self.f = open(filename, 'r')
      # Test the read
      line = self.f.readline().strip()
      if line != 'ITEM: TIMESTEP':
        raise IOError('Not a valid LAMMPS dump file')
      self.f.readline()
      self.f.readline()
      # And get the nb of atoms
      nbatoms = int(self.f.readline())
      self.f.close()
      # If OK, reopen it
      self.f = open(filename, 'r')

    self.sizeconf = 9 + nbatoms

    # Skip useless confs
    if skip:
      line2skip = skip * self.sizeconf

      for _ in xrange(line2skip):
        if self.f.readline() == '':
          raise ValueError('Skipping too much conf')

      logging.info('I skipped %d lines in the Dump for %d confs', line2skip, skip)

    if self.data:
      # Create atom types via data
      i_at = 0
      for atom_type in self.data.atom_types:
        i_at += 1
        atype = molecule.AtomType(str(i_at), atom_type)
        self.atom_types[str(i_at)] = atype

      # Create dummy atoms from data to try to factorize into atomTypes
      for atom in self.data.atoms:
        atype = self.atom_types[str(atom['atom_type_i'])]
        atom = molecule.Atom(str(atom['i']), atom['i'], 1, atom['x'], atom['y'], atom['z'], \
                                 params={'charge': atom['charge']}, atype=atype)

      for atom_type in self.atom_types.values():
        atom_type.factorize()


  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    cnf = {}

    cnf['filename'] = self.filename
    cnf['i'] = self.i
    cnf['dump'] = self
    self.i += 1

    # Reset all atom_type
    for atom_type in self.atom_types.values():
      atom_type.reset()

    line = self.f.next().strip()
    if line == '':
      raise StopIteration
    if line != 'ITEM: TIMESTEP':
      raise ValueError('Not a valid LAMMPS dump file')
    timestep = int(self.f.next())
    cnf['timestep'] = timestep

    line = self.f.next().strip()
    if line != 'ITEM: NUMBER OF ATOMS':
      raise ValueError('Not a valid LAMMPS dump file')
    nbatoms = int(self.f.next())
    cnf['nbatoms'] = nbatoms

    # Box stuff
    line = self.f.next().strip()
    if not line.startswith('ITEM: BOX BOUNDS'):
      raise ValueError('Not a valid LAMMPS dump file')
    cnf['boundaries'] = line.split()[-3:]
    if 'xy xz yz' in line:
      is_tilt = True
    else:
      is_tilt = False

    # The format can't change, so store in raw style
    cnf['box'] = []
    line = self.f.next()
    cnf['box'].append([float(line.split()[0]), float(line.split()[1])])
    if is_tilt:
      xy = float(line.split()[2])

    line = self.f.next()
    cnf['box'].append([float(line.split()[0]), float(line.split()[1])])
    if is_tilt:
      xz = float(line.split()[2])

    line = self.f.next()
    cnf['box'].append([float(line.split()[0]), float(line.split()[1])])
    if is_tilt:
      yz = float(line.split()[2])
    if is_tilt:
      cnf['tilt'] = [xy, xz, yz]
    else:
      cnf['tilt'] = [0.0, 0.0, 0.0]

    # Get item names
    line = self.f.next().strip()
    if not line.startswith('ITEM: ATOMS'):
      raise ValueError('Not a valid LAMMPS dump file')
    fields = line.split()[2:]

    cnf['numpy'] = self.numpy

    s = itertools.islice(self.f, cnf['nbatoms'])
    raw = (list(s), fields, self.sort)
    if self.raw:
      cnf['raw'] = raw
      cnf['atoms'] = None
    else:
      cnf['atoms'] = Dump.raw2atoms(raw, self.numpy)

    return cnf

  @staticmethod
  def raw2atoms(raw, numpy):
    '''Transform raw lines to atoms (array numpy or dict)'''

    (s, fields, sort) = raw
    if numpy:
      #atoms = np.genfromtxt(s, dtype=None, names=fields, usemask=False)
      atoms = pandas.read_csv(StringIO.StringIO(''.join(s)), ' ', header=None, index_col=False, names=fields, as_recarray=True).view(np.ndarray).copy()
      if sort:
        atoms.sort(order='id')
    else:
      atoms = []
      for line in s:
        atom = {}
        values = line.split()
        for (field, value) in zip(fields, values):
          try:
            # Try for int of float
            try:
              atom[field] = int(value)
            except ValueError:
              atom[field] = float(value)
          except ValueError:
            atom[field] = value
        atoms.append(atom)

      if sort:
        atoms = sorted(atoms, key=lambda atom: atom['id'])

    return atoms


class Ave(object):
  '''Class representing a Average file reader'''

  def __init__(self, filename):

    self.filename = filename

    # Try first gzip and fall back to normal if failed
    try:
      self.f = gzip.open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = gzip.open(filename, 'r')
    except IOError:
      self.f = open(filename, 'r')
      # Test the read
      self.f.readline()
      self.f.close()
      # If OK, reopen it
      self.f = open(filename, 'r')

    # Read the first line header
    line = self.f.readline()
    if not line.startswith('#'):
      raise ValueError('File %s does not seem to be a valid LAMMPS ave file' % self.filename)
    self.header = line.split('#', 1)[1].strip()

    # Skip the second useless lines
    line = self.f.readline()
    if not line.startswith('#'):
      raise ValueError('File %s does not seem to be a valid LAMMPS ave file' % self.filename)

    # The third ones contains field names
    line = self.f.readline()
    if not line.startswith('#'):
      raise ValueError('File %s does not seem to be a valid LAMMPS ave file' % self.filename)

    self.fields = line.split('#')[1].split()

  def __iter__(self):
    return self

  def next(self):
    '''The iterator method, returns a dict conf containing all infos of the
       current configuration'''

    cnf = {}

    line = self.f.readline()
    if line == '':
      raise StopIteration

    fields = line.split()
    cnf['timestep'] = int(fields[0])
    nb_line = int(fields[1])

    for field in self.fields:
      cnf[field] = []

    for _ in xrange(nb_line):
      line = self.f.readline()
      fields = line.split()

      # check the coherence between header and line
      if len(fields) != len(self.fields):
        raise ValueError('Invalid data or field names : %s' % self.fields)

      i = 0
      for field in self.fields:
        try:
          val = int(fields[i])
        except ValueError:
          val = float(fields[i])
          val = float(fields[i])
        cnf[field].append(val)
        i += 1

    return cnf

  def to_array(self):
    '''Convert the average to a numpy array'''

    nb_fields = len(self.fields)
    data = list(self)
    types = [('timestep', float)] + zip(self.fields, [float]*nb_fields, [(len(data[0][self.fields[0]]),)] * nb_fields)
    return np.array([(line['timestep'],) + tuple(line[self.fields[i]] for i in xrange(nb_fields)) for line in data], dtype=types)
