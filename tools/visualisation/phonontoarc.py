#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.join(sys.path[0],'../../'))
from fdtoolbox.calculation_set import *
from fdtoolbox.linear_expansion import *
from fdtoolbox.utility import *


try:
  directory = sys.argv[1]
except:
  directory = "."

try:
  mode_num = int(sys.argv[2])
except:
  mode_num = -1

try:
  species = argvtospecies( sys.argv[3] )
except:
  species = None

try:
  usecache = sys.argv[4] != 'nocache'
except:
  usecache = True
  
  
  

  
print '#Reading from directory: %s' % directory
print '#Using cache: %s' % `usecache`

saxis = [0., 0., 1.]
cs=calculation_set.read_directory(directory, saxis, usecache)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

cs.try_fix_displacements()  

lin_exp = linear_expansion( cs )
lin_exp.calculate_expansion_coefficients()

eigstruct = linalg.eig(lin_exp.B_m_n)

if mode_num < 0:
  zeromodes = nonzero(abs(eigstruct[0]) < cs.TRANSLATIONAL_MODE_THRESHOLD)
  eigstruct[0][zeromodes] = 1e9
  mode_num = nonzero(eigstruct[0]==eigstruct[0].min())[0][0]

eigmode = (linalg.eig(lin_exp.B_m_n)[1])[:,mode_num]
print '#Eigenmode number %d with eigenvalue %f'%(mode_num,cs.groundstate.volume*eigstruct[0][mode_num])
sys.stdout.flush()

if species is None:
  species = cs.groundstate.num_atoms*["C"]

cs.groundstate.save_to_arc("/dev/stdout", species)
for i in range(30):
  cs.groundstate.atoms += cos(2*pi*i/30.)/5.*eigmode.reshape((cs.groundstate.num_atoms,3)) 
  cs.groundstate.save_to_arc("/dev/stdout", species, False)
  