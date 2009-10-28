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
  usecache = argv[2] != 'nocache'
except:
  usecache = True
  
print 'Reading from directory: %s' % directory
print 'Using cache: %s' % `usecache`

saxis = [0., 0., 1.]
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

print 'Atom projected charges'
pc = cs.groundstate.projected_charges()
if pc is None:
  print "Not available"
else:
  print mat2str( pc )

bc, units = lin_exp.born_charges()
print 'Born charges (%s)'%units
print mat2str( bc )

ebar, units = lin_exp.piezoelectric_stress_tensor(ionic = False)
print 'Clamped-ion piezoelectric stress tensor in Voigt notation (%s)'%units
print mat2str( tensor2voigt( ebar.T, [1] ) )

etens, units = lin_exp.piezoelectric_stress_tensor()
print 'Full piezoelectric stress tensor in Voigt notation (%s)'%units
print mat2str( tensor2voigt( etens.T, [1] ) )
    
dtens, units = lin_exp.piezoelectric_strain_tensor(ionic = False)
print 'Clamped-ion piezoelectric strain tensor in Voigt notation (%s)'%units
print mat2str( tensor2voigt( dtens.T, [1] ) )

dtens, units = lin_exp.piezoelectric_strain_tensor()
print 'Full piezoelectric strain tensor in Voigt notation (%s)'%units
print mat2str( tensor2voigt( dtens.T, [1] ) )

chitens, units = lin_exp.electric_susceptibility(lattice = False)
print 'Dielectric tensor ionic contribution (%s)'%units
print mat2str( chitens )

chitens, units = lin_exp.electric_susceptibility()
print 'Dielectric tensor ionic and cell contribution (%s)'%units
print mat2str( chitens )
