#!/usr/bin/env python

import sys
import os
from numpy import *

sys.path.append(os.path.join(sys.path[0],'../'))
from fdtoolbox.calculation_set import *


trans_dict = {'RtoC':mat([[1, -1, 1],
                          [-1, 1, 1],
                          [1, 1, -1]]),
              'RtoM':mat([[0.5, 0.5, 0],
                          [0.5,-0.5, 0],
                          [0, 0, 1]])*\
                     mat([[1, -1, 1],
                          [-1, 1, 1],
                          [1, 1, -1]])
              }
inputfile = sys.argv[1]
transformation = sys.argv[2]

if transformation in trans_dict:
  transformation = trans_dict[transformation]
else:
  transformation = mat(transformation)

outputfile = '/dev/stdout'

cell = calculation()
cell.load_from_poscar(inputfile)

newcell = copy.deepcopy(cell)
newcell.unit_cell = transformation*newcell.unit_cell

shape_max = vstack( [ array( abs(transformation).sum(0)), array(transformation.max(0)), [0,0,0] ]).max(0)
shape_min = vstack( [ array(-abs(transformation).sum(0)), array(transformation.min(0)), [0,0,0] ]).min(0)

newcell.atoms = []
newcell.species = []
cur_s = 0
for s_count in cell.species:
  s_count = int(s_count)
  appended_s = 0
  for atom in range( cur_s, cur_s + s_count):
    for i in iterate_all_indices(shape_max, shape_min):
      shifted_atom = cell.atoms[atom,:] + dot(i, cell.unit_cell)
      fract_coords = dot( shifted_atom, newcell.recip_cell )
      if fract_coords.max(1) < 1.0 and fract_coords.min(1) >= 0.0:
        newcell.atoms.extend(array(shifted_atom))
        appended_s += 1

  cur_s += s_count
  newcell.species.append( str(appended_s) )

newcell.atoms = mat(newcell.atoms)


print '# transformed cell a = %.6f, b = %.6f, c = %.6f'%( newcell.cell_a, newcell.cell_b, newcell.cell_c )
print '# alpha = %.4f, beta = %.4f, gamma = %.4f'%( newcell.cell_alpha, newcell.cell_beta, newcell.cell_gamma )

newcell.save_to_poscar('/dev/stdout')
