#!/usr/bin/env python

from fdtoolbox.calculation_set import *
from fdtoolbox.utility import *
from fdtoolbox.linear_expansion import *
import sys
import optparse

parser = optparse.OptionParser(usage='Usage: %prog [options] [DIRECTORY] [EASY_AXIS]')

add_common_options(parser)
parser.add_option("-d", "--decompose", action='store_true', dest='decomp', default=False, help="Show per-mode decomposition")
(options,arguments) = parser.parse_args()

try:
  directory = arguments[0]
except:
  directory = "."
  
try:
  saxis = [float(v) for v in arguments[1].split('/')]  
except:
  saxis = [0., 0., 1.]

if not options.nocache:
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

set_common_options(cs, lin_exp, options)

lin_exp.calculate_expansion_coefficients()


chitens, units = lin_exp.electric_susceptibility(lattice = False)
print 'Dielectric tensor ionic contribution (%s)'%units
print mat2str( chitens )

chitens, units = lin_exp.electric_susceptibility()
print 'Dielectric tensor ionic and cell contribution (%s)'%units
print mat2str( chitens )

if options.decomp:
#  evals = linalg.eig(lin_exp.B_m_n)[0]
  print 'enum\teval\tcontribution'
  for i in range(lin_exp.B_m_n.shape[0]):
    inv_B_m_n, _, ee = safe_inv( lin_exp.B_m_n, cs.TRANSLATIONAL_MODE_THRESHOLD, iterate_components = [i] )
    print '%d\t%f\t%s'%(i, ee[i]*cs.groundstate.volume, mat2str( ( EEAEV_TO_EPSILON0*dot(lin_exp.B_m_alpha.T, dot(inv_B_m_n, lin_exp.B_m_alpha) )).flatten() )[:-1] )

