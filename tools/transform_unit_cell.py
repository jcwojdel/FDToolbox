#!/usr/bin/env python

import sys
from numpy import *
from metoolbox.calculation_set import *


inputfile = sys.argv[1]

outputfile = '/dev/stdout'

cell = calculation()
cell.load_from_poscar(inputfile)

transformation = mat([[1, 1, 0],
                      [1,-1, 0],
                      [0, 0, 1]])*\
                 mat([[1, -1, 1],
                      [-1, 1, 1],
                      [1, 1, -1]])

print cell.cell_a, cell.cell_b, cell.cell_c, cell.cell_alpha, cell.cell_beta, cell.cell_gamma

newcell = copy.deepcopy(cell)
newcell.unit_cell = transformation*newcell.unit_cell

print transformation

print newcell.cell_a, newcell.cell_b, newcell.cell_c, newcell.cell_alpha, newcell.cell_beta, newcell.cell_gamma

