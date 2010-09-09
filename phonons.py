#!/usr/bin/env python

from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *
import sys
import getopt

    
try:                                
  opts, args = getopt.getopt(sys.argv[1:], "zapmc:s", ["zmass","asr","polarity","mpolarity","condition=","sort"])
except getopt.GetoptError:
  print "Usage:"
  print sys.argv[0] + "[options] directory_name easy_axis"
  print "Options:"
  print "-z, --zmass"
  print "\tUse ZMASS from the output in order to calculate real phonons, not just force constant"
  print "\tmatrix eigenvalues."
  print "-a, --asr"
  print "\tForce acoustic sum rule on the FC matrix."
  print "-p, --polarity"
  print "\tPrint out the mode polarities"
  print "-m, --mpolarity"
  print "\tPrint out magnetic strengths of the modes"
  print "-c expr, --condition=expr"
  print "\tPrint out only eigenmodes satisfying given condition."
  print "\tThe condition is a Python expression, which can contain the following variables:"
  print "\t\teigenvalue - mode eigenvalue"
  print "\t\ti - unsorted eigenmode number"
  print "-s, --sort"
  print "\tSort the modes according to their eigenvalue. It also forces printing unsorted mode numbers."
  sys.exit(2)                     


try:
  directory = args[0]
except:
  directory = "."

options = []
condition = 'True'
for opt, arg in opts:
  if opt in ('-z','--zmass'):
    options.append('Z')
  if opt in ('-a','--asr'):
    options.append('A')
  if opt in ('-p','--polarity'):
    options.append('P')
  if opt in ('-m','--mpolarity'):
    options.append('M')
  if opt in ('-s','--sort'):
    options.append('S')
  if opt in ('-c','--condition'):
    condition = arg
    
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

if 'A' in options:
  print 'Enforcing ASR'
  for i in range(0,fc_matrix.shape[0],3):
    asr = zeros([3,3])
    for j in range(0,fc_matrix.shape[0],3):      
      asr += fc_matrix[i:i+3,j:j+3]
      
    fc_matrix[i:i+3,i:i+3] -= asr
    

if 'Z' in options:
  print 'Using atomic masses in the calculations'
  units += '*amu**-1'
  for i in range(len(cs.groundstate.pomass)):
    fc_matrix[i*3:i*3+3,:] /= sqrt(cs.groundstate.pomass[i])
    fc_matrix[:,i*3:i*3+3] /= sqrt(cs.groundstate.pomass[i])
    

evalues, evectors = linalg.eig(fc_matrix)


new_evv = []
for i,eigenvalue in enumerate(evalues):
  if eval(condition):
    new_evv.append((eigenvalue,evectors[:,i],i))
    
if 'S' in options:
  new_evv = sorted(new_evv)
evalues = array([e[0] for e in new_evv])
evectors = hstack([e[1] for e in new_evv])
enumbers = array([e[2] for e in new_evv])

if 'S' in options:
  print 'Unsorted mode numbers:'
  print mat2str(enumbers,'%10d')

if 'Z' in options:
  print 'Phonon frequencies (cm**-1)'
  evalues[:]=[sign(a)*sqrt(abs(a)) for a in evalues]
  print mat2str(SQRTEVANGAMU_TO_CM*evalues.T, '%10.4f')
else:  
  print 'Force constant matrix eigenvalues (%s)'%units
  print mat2str( evalues.T, '%10.4f' )
  
print 'Force constant matrix eigenmodes'
print mat2str( evectors, '%10.4f' )
  
if 'P' in options:   
  print 'Dielectric polarity of the modes'
  print mat2str( dot( lin_exp.born_charges()[0].T, evectors ), '%10.4f' )

if 'M' in options:
  print 'Magnetic "polarity" of the modes (1e-3)'
  print mat2str( 1e3*dot( lin_exp.B_m_mu.T, evectors ), '%10.4f' )
 
