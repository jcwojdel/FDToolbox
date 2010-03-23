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
  polarities = sys.argv[3]
except:
  polarities = ""
  
try:
  saxis = [float(v) for v in sys.argv[4].split('/')]
except:
  saxis = [0., 0., 1.]
    
    
print 'Reading from directory: %s' % directory


cs=calculation_set.read_directory(directory, saxis, True)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')
 
calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5
cs.try_fix_displacements()
cs.try_fix_polarizations()

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
  
if 'P' in polarities:   
  print 'Dielectric polarity of the modes'
  print mat2str( dot( lin_exp.born_charges()[0].T, evectors ), '%10.4f' )

if 'M' in polarities:
  print 'Magnetic "polarity" of the modes (1e-3)'
  print mat2str( 1e3*dot( lin_exp.B_m_mu.T, evectors ), '%10.4f' )
 

#print 'Force constant matrix eigenvalues (%s)'%units
#print mat2str( evalues.T, '%10.4f' )
#print 'Force constant matrix eigenmodes'
#print mat2str( evectors, '%10.4f' )
