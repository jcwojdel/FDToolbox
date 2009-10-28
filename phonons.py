#!/usr/bin/env python

from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *
import sys

try:
  directory = sys.argv[1]
except:
  directory = "."

try:
  condition = sys.argv[2]
except:
  condition = "True"
  
try:
  usecache = argv[3] != 'nocache'
except:
  usecache = True
  
print 'Reading from directory: %s' % directory
print 'Using cache: %s' % `usecache`

saxis = [0., 0., 1.]
cs=calculation_set.read_directory(directory, saxis, usecache)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')
 
cs.try_fix_displacements()

lin_exp = linear_expansion( cs )
lin_exp.calculate_expansion_coefficients()

fc_matrix, units = lin_exp.force_constant_matrix()
evalues, evectors = linalg.eig(fc_matrix)


new_evalues=[]
new_evectors=[]
for i,eigenvalue in enumerate(evalues):
  if eval(condition):
    new_evalues.append(eigenvalue)
    new_evectors.append(evectors[:,i])
    
evalues = array(new_evalues)
evectors = hstack(new_evectors)

print 'Force constant matrix eigenvalues (%s)'%units
print mat2str( evalues.T, '%10.4f' )
print 'Force constant matrix eigenmodes'
print mat2str( evectors, '%10.4f' )
    


#print 'Force constant matrix eigenvalues (%s)'%units
#print mat2str( evalues.T, '%10.4f' )
#print 'Force constant matrix eigenmodes'
#print mat2str( evectors, '%10.4f' )
