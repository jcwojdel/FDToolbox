#!/usr/bin/env python

from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *
import sys
import optparse

parser = optparse.OptionParser(usage='Usage: %prog [options] [DIRECTORY] [EASY_AXIS]')

parser.add_option("-c", "--condition", action='store', type='string', dest='condition', help="Print out only eigenmodes satisfying given condition. The condition EXPR must be a Python expression, which can contain the following variables: eigenvalue - mode eigenvalue, and i - unsorted eigenmode number.", metavar="EXPR")
parser.add_option("-z", "--zmass", action='store_true', dest='zmass', default=False, help="Use ZMASS from the output in order to calculate real phonons, not just force constant matrix eigenvalues.")
parser.add_option("-a", "--asr", action='store_true', dest='asr', default=False, help="Force acoustic sum rule on the FC matrix.")
parser.add_option("-p", "--polarity", action='store_true', dest='pol', default=False, help="Print out the mode polarities.")
parser.add_option("-m", "--mpolarity", action='store_true', dest='mpol', default=False, help="Print out magnetic strengths of the modes.")
parser.add_option("-s", "--sort", action='store_true', dest='sort', default=False, help="Sort the modes according to their eigenvalue. It also forces printing unsorted mode numbers.")
parser.add_option("-v", "--vectors", action='store_true', dest='vec', default=False, help="Print out the eigen vectors, not only their eigenvalues.")

(options,arguments) = parser.parse_args()

if options.condition is None:
  options.condition = 'True'
    
try:
  directory = arguments[0]
except:
  directory = "."
    
try:
  saxis = [float(v) for v in args[1].split('/')]
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

if options.asr:
  print 'Enforcing ASR'
  for i in range(0,fc_matrix.shape[0],3):
    asr = zeros([3,3])
    for j in range(0,fc_matrix.shape[0],3):      
      asr += fc_matrix[i:i+3,j:j+3]
      
    fc_matrix[i:i+3,i:i+3] -= asr
    

if options.zmass:
  print 'Using atomic masses in the calculations'
  units += '*amu**-1'
  for i in range(len(cs.groundstate.pomass)):
    fc_matrix[i*3:i*3+3,:] /= sqrt(cs.groundstate.pomass[i])
    fc_matrix[:,i*3:i*3+3] /= sqrt(cs.groundstate.pomass[i])
    

evalues, evectors = linalg.eig(fc_matrix)


new_evv = []

if options.zmass:
  evalues[:]=[SQRTEVANGAMU_TO_CM*sign(a)*sqrt(abs(a)) for a in evalues]
  
for i,eigenvalue in enumerate(evalues):
  if eval(options.condition):
    new_evv.append((eigenvalue,evectors[:,i],i))
    
if options.sort:
  new_evv = sorted(new_evv)
evalues = array([e[0] for e in new_evv])
evectors = hstack([e[1] for e in new_evv])
enumbers = array([e[2] for e in new_evv])

if options.sort:
  print 'Unsorted mode numbers:'
  print mat2str(enumbers,'%10d')

if options.zmass:
  print 'Phonon frequencies (cm**-1)'  
  print mat2str(evalues.T, '%10.4f')
else:  
  print 'Force constant matrix eigenvalues (%s)'%units
  print mat2str( evalues.T, '%10.4f' )

if options.vec:  
  print 'Force constant matrix eigenmodes'
  print mat2str( evectors, '%10.4f' )
  
if options.pol:   
  print 'Dielectric polarity of the modes'
  print mat2str( dot( lin_exp.born_charges()[0].T, evectors ), '%10.4f' )

if options.mpol:
  print 'Magnetic "polarity" of the modes (1e-3)'
  print mat2str( 1e3*dot( lin_exp.B_m_mu.T, evectors ), '%10.4f' )
 
