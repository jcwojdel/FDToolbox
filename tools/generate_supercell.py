#!/usr/bin/env python

import sys
from numpy import *
from metoolbox.calculation_set import *

if len(sys.argv) < 3:
  name = (sys.argv[0].split('/'))[-1]
  print 'Usage:\n\
  %s inputfilename repetitions {outputfilename}\n\
  filename - name of the POSCAR/CONTCAR file containing optimised system\n\
  repetitions - number of cell repetitions in each direction\n\
Example:\n\
  %s ./POSCAR 2x2x2 \n' % ( name, name )

inputfile = sys.argv[1]
repetitions = [int(i) for i in sys.argv[2].split('x')]

if len(sys.argv) > 3:
  outputfile = sys.argv[3]
else:
  outputfile = '/dev/stdout'

cell = calculation()
cell.load_from_poscar(inputfile)

cell.supercell(repetitions).save_to_poscar(outputfile)

