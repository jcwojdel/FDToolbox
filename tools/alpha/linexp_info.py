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
  directory = sys.argv[1]
except:
  directory = "."
  
try:
  saxis = [float(v) for v in sys.argv[2].split('/')]  
except:
  saxis = [0., 0., 1.]

usecache = True
  
print 'Reading from directory: %s' % directory
print 'Using saxis: %s' % `saxis`
print 'Using cache: %s' % `usecache`

cs=calculation_set.read_directory(directory, saxis, usecache)

cs.set_ionic('.*/calc_.._.*[0123]')
cs.set_lattice('.*/calc_00_...')
cs.set_groundstate('.*/calc_00_0')

calculation_set.BPH_MAXDIST = 0.25
calculation_set.BPH_CORRECTION = 0.5

cs.try_fix_displacements()  
cs.try_fix_polarizations()

lin_exp = linear_expansion( cs )
lin_exp.calculate_expansion_coefficients()

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
