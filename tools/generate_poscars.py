#!/usr/bin/env python

from numpy import *


from metoolbox.calculation_set import calculation
import copy
import sys
import os

disp_length = [0.005, 0.005, 0.005]
lat_disp_length = [0.01, 0.01, 0.01]

def move_all_atoms(calc):
  yield 'calc_00_0', calc
  for numat, at in enumerate(calc.atoms):
    for direction in range(3):
      displacement = zeros((1,3))
      displacement[0, direction] = disp_length[direction % len(disp_length)]
      displacement = displacement*calc.recip_cell
      
      for orientation in ['+', '-']:
        c = copy.deepcopy(calc)
        c.atoms[numat] += displacement
        yield 'calc_%02d_%s%d' % (numat+1, orientation, direction+1), c
      
        displacement *= -1

def move_all_lattice(calc):
  for axnum, axname in enumerate(['a', 'b', 'c']):
    for dirnum, dirname in enumerate(['X', 'Y', 'Z']):
      displacement = zeros((1,3))
      displacement[0, dirnum] = lat_disp_length[dirnum % len(lat_disp_length)]
      for orientation in ['+', '-']:
        c = copy.deepcopy(calc)
        c.unit_cell[axnum] += displacement
        c.atoms = dot( calc.atoms, dot( linalg.inv( calc.unit_cell ), c.unit_cell ) )
        yield 'calc_00_%s%s%s' % (axname, orientation, dirname), c
      
        displacement *= -1

def generate_all_calcs(calc):
  for name, c in move_all_atoms(calc):
    yield name, c
  for name, c in move_all_lattice(calc):
    yield name, c
  
if len(sys.argv) < 2:
  name = (sys.argv[0].split('/'))[-1]
  print 'Usage:\n\
  %s filename {ion_displacement_lengths {lattice_displacement_lengths}}\n\
  filename - name of the POSCAR/CONTCAR file containing optimised system\n\
  ion_displacement_lengths - ion displacement in angstroms single value for uniform displacement,\n\
                             or quoted triplet of values\n\
  lattice_displacement_lengths - same as ion_displacement_lengths but for lattice vectors\n\
Example:\n\
  %s ./POSCAR 0.005 "0.01 0.01 0.05"\n' % ( name, name )
  sys.exit()
      
system = calculation()
system.load_from_poscar(sys.argv[1])


if len(sys.argv) > 2:
  disp_length = [ float(d) for d in sys.argv[2].split() ]

if len(sys.argv) > 3:
  lat_disp_length = [ float(d) for d in sys.argv[3].split() ]

generator = generate_all_calcs
if max(disp_length) == 0.:
  generator = move_all_lattice
if max(lat_disp_length) == 0.:
  generator = move_all_atoms

for name, calc in generator(system):
  try:
    os.mkdir(name)
  except:
    pass
  

  calc.save_to_poscar('/'.join( [name, 'POSCAR'] ))
  print name
  