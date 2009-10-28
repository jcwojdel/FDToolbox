from metoolbox.calculation_set import *
from metoolbox.utility import *
from metoolbox.linear_expansion import *
import sys
import os

try:
  path = sys.argv[1]
except:
  path = "."

try:
  calculation.saxis = [float(v) for v in sys.argv[-1].split('/')]  
except:
  calculation.saxis = [0., 0., 1.]

cs = calculation_set()
  
for path in sys.argv[1:]:  
  if os.path.isfile(path) or os.path.isfile(path+'.gz'):
    print 'Reading from single file: %s' % path
    calc=calculation()
    cs.add_calculation( calc )
  elif os.path.isdir(path):
    if os.path.isfile( os.path.join(path, 'OUTCAR') ) or os.path.isfile( os.path.join(path, 'OUTCAR.gz') ):
      print 'Reading single calculation from directory: %s' % path
      calc=calculation()
      calc.load_from_outcar( os.path.join(path, 'OUTCAR') )
      cs.add_calculation( calc )      
    else:
      print 'Reading calculations from directory: %s' % path 
      cs = calculation_set.read_directory(path, calculation.saxis, True)
  else:
    pass

cs.set_ionic('.')
cs.groundstate = cs.calculations()[0]
cs.try_polarizations_chain()  
#cs.try_fix_polarizations()

for calc in cs.calculations():
  print 'Calc name: %s' % calc.name
  print 'Calc fileID: %s' % calc.fileID
  print 'magnetization: %s' % mat2str(calc.magnetization)
  print 'polarization: %s' % mat2str(calc.polarization)
