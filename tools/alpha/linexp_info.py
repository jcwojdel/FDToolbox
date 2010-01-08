#!/usr/bin/env python

import sys
import os

sys.path.append(os.path.join(sys.path[0],'../../'))
from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *

loggable.log_level = LOG_ALLINFO
calculation_set.TRANSLATIONAL_MODE_THRESHOLD = 0.001

try:
  directory1 = sys.argv[1]
except:
  directory1 = "."

try:
  directory2 = sys.argv[2]
except:
  directory2 = directory1

try:
  fraction = float(sys.argv[3])
except:
  fraction = 0.5

try:
  saxis = [float(v) for v in sys.argv[4].split('/')]  
except:
  saxis = [0., 0., 1.]
  
  

usecache = True
  
print 'Reading from directories: %s %s' % (directory1, directory2)
print 'Interpolation: %f' % fraction
print 'Using saxis: %s' % `saxis`
print 'Using cache: %s' % `usecache`

cs1=calculation_set.read_directory(directory1, saxis, usecache)

cs1.set_ionic('.*/calc_.._.*[0123]')
cs1.set_lattice('.*/calc_00_...')
cs1.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs1.try_fix_displacements()  
cs1.try_fix_polarizations()

cs2=calculation_set.read_directory(directory2, saxis, usecache)

cs2.set_ionic('.*/calc_.._.*[0123]')
cs2.set_lattice('.*/calc_00_...')
cs2.set_groundstate('.*/calc_00_0')

cs2.try_fix_displacements()  
cs2.try_fix_polarizations()

cs1.add_lattice_constraint( array([[1,0,0],[0,1,0]]), 'xyshear' )
cs2.add_lattice_constraint( array([[1,0,0],[0,1,0]]), 'xyshear' )
lin_exp1 = linear_expansion( cs1 )
lin_exp1.calculate_expansion_coefficients()

lin_exp2 = linear_expansion( cs2 )
lin_exp2.calculate_expansion_coefficients()

def interpolate(a, b, f):
  return a*(1.-f)+b*f

lin_exp = linear_expansion( cs1 )
for quantity in [ "volume", "B_alpha_beta", "B_alpha_mu", "B_mu_nu", "B_m_n", "B_j_k", "B_m_j", "B_m_alpha", "B_alpha_j", "B_m_mu", "B_mu_j" ]:
  lin_exp.__dict__[quantity] = interpolate( lin_exp1.__dict__[quantity], lin_exp2.__dict__[quantity], fraction)
  

lin_exp.calculate_expansion_coefficients(update_from_calcset=False)


"""
print "B_m_n"
print mat2str(lin_exp.B_m_n, "%8E")
  
print "B_j_k"
print mat2str(lin_exp.B_j_k, "%8E")

print "B_m_j"
print mat2str(lin_exp.B_m_j, "%8E")
  
print "B_m_alpha"
print mat2str(lin_exp.B_m_alpha, "%8E")

print "B_alpha_j"
print mat2str(lin_exp.B_alpha_j, "%8E")

print "B_m_mu"
print mat2str(lin_exp.B_m_mu, "%8E")

print "B_mu_j"
print mat2str(lin_exp.B_mu_j, "%8E")
"""

#print mat2str(lin_exp.magneto_electric_coupling()[0].T.flatten())
#print mat2str( linalg.eig(lin_exp.elastic_tensor()[0]) [0] ) 
#print mat2str(1e8*lin_exp.piezomagnetic_strain_tensor()[0].T.flatten())
print mat2str( linalg.eig(lin_exp.elastic_tensor(ionic=False)[0]) [0] )