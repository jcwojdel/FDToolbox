#!/usr/bin/env python

from metoolbox.calculation_set import *
from metoolbox.utility import *
from metoolbox.linear_expansion import *
import sys

try:
  directory = sys.argv[1]
except:
  directory = "."
  
try:
  saxis = [float(v) for v in sys.argv[2].split('/')]  
except:
  saxis = [0., 0., 1.]

try:
  usecache = argv[3] != 'nocache'
except:
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


print 'Total magnetization of groundstate'
print mat2str( cs.groundstate.magnetization )

print 'Atom projected magnetizations'
pm = cs.groundstate.projected_magnetizations()
if pm is None:
  print "Not available"
else:
  print mat2str( pm )

ebar, units = lin_exp.piezomagnetic_stress_tensor(ionic = False)
print 'Clamped-ion piezomagnetic stress tensor in Voigt notation (1e4 %s)'%units
print mat2str( 1e-4*tensor2voigt( ebar.T, [1] ) )

etens, units = lin_exp.piezomagnetic_stress_tensor()
print 'Full piezomagnetic stress tensor in Voigt notation (1e4 %s)'%units
print mat2str( 1e-4*tensor2voigt( etens.T, [1] ) )
    
dtens, units = lin_exp.piezomagnetic_strain_tensor(ionic = False)
print 'Clamped-ion piezomagnetic strain tensor in Voigt notation (1e-8 %s)'%units
print mat2str( 1e8*tensor2voigt( dtens.T, [1] ) )

dtens, units = lin_exp.piezomagnetic_strain_tensor()
print 'Full piezomagnetic strain tensor in Voigt notation (1e-8 %s)'%units
print mat2str( 1e8*tensor2voigt( dtens.T, [1] ) )

chitens, units = lin_exp.magnetic_susceptibility(lattice = False)
print 'Magnetic susceptibility tensor ionic contribution (%s)'%units
print mat2str( chitens )

chitens, units = lin_exp.magnetic_susceptibility()
print 'Magnetic susceptibility tensor ionic and cell contribution (%s)'%units
print mat2str( chitens )

