#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.join(sys.path[0],'../../'))
from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *

species = argvtospecies( sys.argv[1] )

cs = calculation_set()

loggable.log_level = LOG_ERROR
  
for path in sys.argv[2:]:  
  if os.path.isfile(path) or os.path.isfile(path+'.gz'):
    #print 'Reading from single file: %s' % path
    calc=calculation()
    cs.add_calculation( calc )
  elif os.path.isdir(path):
    if os.path.isfile( os.path.join(path, 'OUTCAR') ) or os.path.isfile( os.path.join(path, 'OUTCAR.gz') ):
      #print 'Reading single calculation from directory: %s' % path
      calc=calculation()
      calc.load_from_outcar( os.path.join(path, 'OUTCAR') )
      cs.add_calculation( calc )      
    else:
      #print 'Reading calculations from directory: %s' % path 
      cs = calculation_set.read_directory(path, calculation.saxis, True)
  else:
    pass

cs.set_ionic('.')
cs.groundstate = cs.calculations()[0]

newfile=True
for calc in cs.calculations():
  calc.save_to_arc("/dev/stdout", species, newfile)
  newfile = False
  